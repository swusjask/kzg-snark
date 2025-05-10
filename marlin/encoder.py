import sys, os
from sage.all import vector, prod, PolynomialRing, GF

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fft_ff import fft_ff_interpolation


class Encoder:
    """
    Handles the encoding of R1CS constraint systems into polynomial form for Marlin protocol.
    
    This class provides methods to convert R1CS matrices, witness vectors, and linear combinations
    into their polynomial representations as required by the Marlin protocol.
    """
    
    def __init__(self, q):
        """
        Initialize the R1CS encoder with the specified field.

        Args:
            q: Prime field size
        """
        self.Fq = GF(q)  # Finite field GF(curve_order)
        self.R = PolynomialRing(self.Fq, 'X')
        self.X = self.R.gen()

    def update_state(self, A, B, C):
        """
        Update the encoder state with R1CS constraint matrices.

        Args:
            A, B, C: The R1CS constraint matrices
        """
        self.A = A
        self.B = B
        self.C = C
        
        # Calculate appropriate subgroup size for domain H (must be power of 2)
        self.n = self.find_subgroup_size(max(A.nrows(), A.ncols()))

        # Calculate appropriate subgroup size for domain K (must be power of 2)
        num_nonzero_A = len(A.nonzero_positions())
        num_nonzero_B = len(B.nonzero_positions())
        num_nonzero_C = len(C.nonzero_positions())
        self.m = self.find_subgroup_size(
            max(num_nonzero_A, num_nonzero_B, num_nonzero_C)
        )

        # Generate multiplicative subgroups H and K
        self.g_H = self.Fq(1).nth_root(self.n)
        self.g_K = self.Fq(1).nth_root(self.m)
        self.H = [self.g_H**i for i in range(self.n)]
        self.K = [self.g_K**i for i in range(self.m)]
        
        # Compute vanishing polynomials
        self.v_H = self.X**self.n - 1
        self.v_K = self.X**self.m - 1

    def find_subgroup_size(self, n):
        """
        Find the smallest power of 2 that's greater than or equal to n.

        Args:
            n: The minimum required size

        Returns:
            Power of 2 that is >= n
        """
        return 2 ** ((n - 1).bit_length())

    def u_H(self, a, b):
        """
        Compute the polynomial u_H(a, b) = (v_H(a) - v_H(b))/(a - b).
        
        This is a key component for implementing the holographic proof optimization.

        Args:
            a, b: Field elements
            
        Returns:
            Value of u_H(a, b)
        """
        if a == b:
            # When a=b, u_H(a,a) is the formal derivative of v_H at a
            return self.v_H.derivative()(a)
        else:
            return (self.v_H(a) - self.v_H(b)) / (a - b)

    def encode_matrices(self):
        """
        Encode the R1CS matrices A, B, C into polynomial form for Marlin.

        This method converts each constraint matrix into row, column, and value polynomials
        following the Marlin protocol specification.
        
        Returns:
            dict: Encoded polynomials for each matrix component
        """
        # Precompute u_H values for efficiency
        u_H_diag = {h: self.u_H(h, h) for h in self.H}

        # Dictionary to store all encoded polynomials
        encoded = {}
        
        # Encode each constraint matrix (A, B, C)
        for name, M in [("A", self.A), ("B", self.B), ("C", self.C)]:
            # Get non-zero positions and values
            nonzero_positions = list(M.nonzero_positions())

            # Prepare arrays for row, column, and values
            row_values = [self.Fq(0)] * self.m
            col_values = [self.Fq(0)] * self.m
            val_values = [self.Fq(0)] * self.m

            # Prepare data for K-interpolation
            for k, (i, j) in enumerate(nonzero_positions):
                row_values[k] = self.H[i]
                col_values[k] = self.H[j]
                # Adjust value by u_H factors as required by Marlin protocol
                val_values[k] = self.Fq(M[i, j]) / (
                    u_H_diag[self.H[i]] * u_H_diag[self.H[j]]
                )

            # Use FFT interpolation to get polynomials
            row_poly = fft_ff_interpolation(row_values, self.g_K, self.Fq)
            col_poly = fft_ff_interpolation(col_values, self.g_K, self.Fq)
            val_poly = fft_ff_interpolation(val_values, self.g_K, self.Fq)

            # Store in dictionary with descriptive keys
            encoded[f"row_{name}"] = row_poly
            encoded[f"col_{name}"] = col_poly
            encoded[f"val_{name}"] = val_poly

        return encoded

    def encode_witness(self, z, x_size):
        """
        Encode the witness vector into polynomials, splitting it into public input and private witness.

        Args:
            z: The full variable assignment vector
            x_size: Size of the public input
            
        Returns:
            dict: Encoded polynomials for the witness
        """
        # Convert z to field elements
        z = [self.Fq(zi) for zi in z]

        # Split z into x (public input) and w (private witness)
        x, w = z[:x_size], z[x_size:]

        X = self.X

        # Create Lagrange polynomial for x (public input)
        x_points = [(self.H[i], x[i]) for i in range(len(x))]
        x_poly = self.R.lagrange_polynomial(x_points)

        # Calculate vanishing polynomial for x positions
        v_H_x = prod([X - self.H[i] for i in range(len(x))])

        # Encode w using the optimization approach
        # First, create values array with zeros for x positions
        values = [self.Fq(0)] * len(x)

        # Add witness values adjusted by public input polynomial
        for i, wi in enumerate(w):
            values.append(wi - x_poly(self.H[i + len(x)]))

        # Pad values to match H size if needed
        padding_size = self.n - len(values)
        if padding_size > 0:
            values.extend([self.Fq(0)] * padding_size)

        # Interpolate to get polynomial f
        f = fft_ff_interpolation(values, self.g_H, self.Fq)

        # Calculate w_poly = f / v_H_x
        w_poly = f // v_H_x
        assert w_poly * v_H_x == f, "w_poly is not well-defined"

        # Reconstruct z_poly = w_poly * v_H_x + x_poly
        z_poly = w_poly * v_H_x + x_poly

        return {
            "x_poly": x_poly,
            "w_poly": w_poly,
            "z_poly": z_poly,
            "x": x,
            "w": w,
        }

    def encode_linear_combinations(self, z):
        """
        Encode the linear combinations zA, zB, zC into polynomial form.

        Args:
            z: The full variable assignment vector
            
        Returns:
            dict: Encoded polynomials for the linear combinations
        """
        # Convert z to a vector of field elements
        z_vector = vector(self.Fq, [self.Fq(zi) for zi in z])

        # Compute the linear combinations
        zA = self.A * z_vector
        zB = self.B * z_vector
        zC = self.C * z_vector

        # Convert to lists
        zA_list = list(zA)
        zB_list = list(zB)
        zC_list = list(zC)

        # Pad lists to match subgroup size
        if len(zA_list) < self.n:
            zA_list.extend([self.Fq(0)] * (self.n - len(zA_list)))
        if len(zB_list) < self.n:
            zB_list.extend([self.Fq(0)] * (self.n - len(zB_list)))
        if len(zC_list) < self.n:
            zC_list.extend([self.Fq(0)] * (self.n - len(zC_list)))

        # Use FFT interpolation to get polynomials
        zA_poly = fft_ff_interpolation(zA_list, self.g_H, self.Fq)
        zB_poly = fft_ff_interpolation(zB_list, self.g_H, self.Fq)
        zC_poly = fft_ff_interpolation(zC_list, self.g_H, self.Fq)

        return {
            "zA_poly": zA_poly,
            "zB_poly": zB_poly,
            "zC_poly": zC_poly,
            "zA": zA_list,
            "zB": zB_list,
            "zC": zC_list,
        }
    
if __name__ == "__main__":
    from py_ecc.optimized_bn128 import curve_order
    import random, pickle

    print("Testing Marlin Encoder")
    print("=" * 60)

    # Load test instance
    with open("R1CS_INSTANCE.pkl", "rb") as f:
        RICS_INSTANCE = pickle.load(f)

    A, B, C, z = RICS_INSTANCE["A"], RICS_INSTANCE["B"], RICS_INSTANCE["C"], RICS_INSTANCE["z"]

    # Initialize encoder
    encoder = Encoder(curve_order)
    encoder.update_state(A, B, C)
    
    # Encode matrices
    print("\nEncoding R1CS matrices...")
    encoded_matrices = encoder.encode_matrices()
    
    # Encode witness (assuming first 5 elements are public inputs)
    print("Encoding witness...")
    x_size = 5
    encoded_witness = encoder.encode_witness(z, x_size)
    
    # Encode linear combinations
    print("Encoding linear combinations...")
    encoded_combinations = encoder.encode_linear_combinations(z)
    
    # Test polynomial relationships on random elements from H
    num_tests = 5
    all_tests_passed = True
    
    print("\nTesting polynomial relationships at random points...")
    for i in range(num_tests):
        # Pick a random element from H
        kappa = random.choice(encoder.H)
        
        # Test: zA(κ) * zB(κ) - zC(κ) = 0
        zA_kappa = encoded_combinations["zA_poly"](kappa)
        zB_kappa = encoded_combinations["zB_poly"](kappa)
        zC_kappa = encoded_combinations["zC_poly"](kappa)
        
        entry_wise_product = (zA_kappa * zB_kappa - zC_kappa)
        test_passed = entry_wise_product == 0
        
        if not test_passed:
            print(f"❌ Test {i+1}: Entry-wise product equality failed at κ={kappa}")
            all_tests_passed = False
            print(f"   zA({kappa}) * zB({kappa}) - zC({kappa}) = {entry_wise_product} ≠ 0")
    
    if all_tests_passed:
        print("✅ All polynomial relationship tests passed!")
    
    print("\nEncoder output statistics:")
    print(f"- Number of encoded matrices: 3 (A, B, C)")
    print(f"- Number of matrix polynomials: {len(encoded_matrices)}")
    print(f"- Subgroup sizes: H: {encoder.n}, K: {encoder.m}")
    
    print("\n✅ Marlin Encoder test completed successfully!")
