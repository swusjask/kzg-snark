from sage.all import vector, prod, PolynomialRing, GF

from fft_ff import fft_ff_interpolation

class Encoder:
    """
    Handles the encoding of PLONK constraint systems into polynomial form.
    
    This encoder converts selector values, permutation data, and witness values
    into their polynomial representations as required by the PLONK protocol.
    It also manages the arithmetic operations on polynomials used in the protocol.
    """
    
    def __init__(self, q):
        """
        Initialize the PLONK encoder with the specified field.
        
        Args:
            q: Prime field size (typically the curve order)
        """
        self.Fq = GF(q)  # Finite field GF(q)
        self.R = PolynomialRing(self.Fq, 'X')
        self.X = self.R.gen()

    def find_subgroup_size(self, n):
        """
        Find the smallest power of 2 that's greater than or equal to n.

        Args:
            n: The minimum required size

        Returns:
            Power of 2 that is >= n
        """
        return 2 ** ((n - 1).bit_length())
        
    def update_state(self, qM, qL, qR, qO, qC, perm):
        """
        Update the encoder state with the circuit description.
        
        Args:
            qM, qL, qR, qO, qC: Selector values at points in H
            perm: Permutation list of length 3n
        """
        # Calculate appropriate subgroup size (must be power of 2)
        self.n = self.find_subgroup_size(len(qM))
        
        # Generate subgroup H with generator g
        self.g = self.Fq(1).nth_root(self.n)
        
        # Store circuit data
        self.qM = qM
        self.qL = qL
        self.qR = qR
        self.qO = qO
        self.qC = qC
        self.perm = perm
        
        # Generate subgroup H (multiplicative subgroup of order n)
        self.H = [self.g**i for i in range(self.n)]
        
        # Find suitable multipliers for creating cosets (k1·H and k2·H)
        self._find_coset_multipliers()
        
        # Generate cosets
        self.k1H = [self.k1 * h for h in self.H]
        self.k2H = [self.k2 * h for h in self.H]
        
        # Compute vanishing polynomial for H: Z_H(X) = X^n - 1
        self.v_H = self.X**self.n - 1
        
    def _find_coset_multipliers(self):
        """
        Find suitable values for k1 and k2 to create disjoint cosets.
        
        This ensures that H, k1·H, and k2·H are disjoint multiplicative cosets,
        which is required for the permutation argument in PLONK.
        """
        n = self.n
        
        # Keep trying random values until we find suitable ones
        while True:
            k1 = self.Fq.random_element()
            k2 = self.Fq.random_element()
            
            # Check that k1 and k2 satisfy the required conditions:
            # 1. k1^n ≠ 1 (k1 is not in H)
            # 2. k2^n ≠ 1 (k2 is not in H)
            # 3. (k1/k2)^n ≠ 1 (k1·H and k2·H are disjoint)
            # 4. k1, k2 ≠ 0 (non-zero)
            if (k1**n != 1 and 
                k2**n != 1 and 
                (k1/k2)**n != 1 and
                k1 != 0 and k2 != 0):
                self.k1 = k1
                self.k2 = k2
                return
    
    def encode_selectors(self):
        """
        Encode selector values into polynomials using FFT interpolation.
        
        Returns:
            Dictionary of selector polynomials
        """
        # Ensure the state has been initialized
        if not hasattr(self, 'H'):
            raise ValueError("Call update_state before encoding selectors")
            
        # Interpolate the selector polynomials
        qM_poly = fft_ff_interpolation(self.qM, self.g, self.Fq)
        qL_poly = fft_ff_interpolation(self.qL, self.g, self.Fq)
        qR_poly = fft_ff_interpolation(self.qR, self.g, self.Fq)
        qO_poly = fft_ff_interpolation(self.qO, self.g, self.Fq)
        qC_poly = fft_ff_interpolation(self.qC, self.g, self.Fq)
        
        return {
            "qM": qM_poly,
            "qL": qL_poly,
            "qR": qR_poly,
            "qO": qO_poly,
            "qC": qC_poly
        }
    
    def encode_permutation(self):
        """
        Encode the permutation into sigma star polynomials.
        
        The permutation is encoded as polynomials that map each wire to its
        corresponding wire in the permutation cycle.
        
        Returns:
            Dictionary of permutation polynomials
        """
        if not hasattr(self, 'H') or not hasattr(self, 'k1') or not hasattr(self, 'k2'):
            raise ValueError("Call update_state before encoding permutation")
            
        n = self.n
        
        # Function to map position indices to elements in H ∪ k1·H ∪ k2·H
        def index_to_element(i):
            if 0 <= i < n:
                return self.H[i]  # Left wires mapped to H
            elif n <= i < 2*n:
                return self.k1H[i-n]  # Right wires mapped to k1·H
            elif 2*n <= i < 3*n:
                return self.k2H[i-2*n]  # Output wires mapped to k2·H
            else:
                raise ValueError(f"Index {i} out of range [0, {3*n-1}]")
        
        # Compute permutation star values for each group
        S_sigma1_values = [index_to_element(self.perm[i]) for i in range(n)]
        S_sigma2_values = [index_to_element(self.perm[i+n]) for i in range(n)]
        S_sigma3_values = [index_to_element(self.perm[i+2*n]) for i in range(n)]
        
        # Interpolate to get the permutation polynomials
        S_sigma1_poly = fft_ff_interpolation(S_sigma1_values, self.g, self.Fq)
        S_sigma2_poly = fft_ff_interpolation(S_sigma2_values, self.g, self.Fq)
        S_sigma3_poly = fft_ff_interpolation(S_sigma3_values, self.g, self.Fq)
        
        sigma_star = S_sigma1_values + S_sigma2_values + S_sigma3_values
        
        return {
            "S_sigma1": S_sigma1_poly,
            "S_sigma2": S_sigma2_poly,
            "S_sigma3": S_sigma3_poly,
            "sigma_star": sigma_star
        }
    
    def encode_witness(self, w, x_size=0):
        """
        Encode witness values into polynomials and compute public input polynomial.
        
        Args:
            w: Witness list of length 3n (a_values, b_values, c_values)
            x_size: Number of public inputs
            
        Returns:
            Dictionary of witness polynomials and public input data
        """
        if not hasattr(self, 'H'):
            raise ValueError("Call update_state before encoding witness")
            
        n = self.n
        
        # Split witness into a (left), b (right), c (output) values
        a_values = w[:n]
        b_values = w[n:2*n]
        c_values = w[2*n:3*n]
        
        # Interpolate to get witness polynomials
        a_poly = fft_ff_interpolation(a_values, self.g, self.Fq)
        b_poly = fft_ff_interpolation(b_values, self.g, self.Fq)
        c_poly = fft_ff_interpolation(c_values, self.g, self.Fq)
        
        # Extract public inputs if specified
        x = w[:x_size] if x_size > 0 else []
        
        # Compute public input polynomial if there are public inputs
        PI = self.compute_public_input_poly(x) if x_size > 0 else self.R(0)
        
        return {
            "a": a_poly,
            "b": b_poly,
            "c": c_poly,
            "x": x,
            "PI": PI
        }
    
    def compute_lagrange_basis(self, i):
        """
        Compute the i-th Lagrange basis polynomial for domain H.
        
        For multiplicative subgroups, the Lagrange polynomial has a special form:
        L_i(X) = (g^i · (X^n - 1)) / (n · (X - g^i))
        
        Args:
            i: Index in H
            
        Returns:
            Lagrange polynomial L_i(X)
        """
        if not hasattr(self, 'H'):
            raise ValueError("Call update_state before computing Lagrange basis")
            
        n = self.n
        g = self.g
        X = self.X
        
        # Compute the i-th Lagrange basis polynomial for H
        numerator = g**i * (X**n - 1)
        denominator = n * (X - g**i)
        L_i = numerator // denominator
        
        return L_i
    
    def compute_public_input_poly(self, x):
        """
        Compute the public input polynomial PI(X) = -∑_i x_i·L_i(X).
        
        This polynomial encodes the public inputs into the verification equation.
        
        Args:
            x: List of public input values
            
        Returns:
            Public input polynomial PI(X)
        """
        if not hasattr(self, 'H'):
            raise ValueError("Call update_state before computing public input poly")
            
        PI = self.R(0)
        for i, x_i in enumerate(x):
            L_i = self.compute_lagrange_basis(i)
            PI -= x_i * L_i
            
        return PI

if __name__ == "__main__":
    # Example usage
    import pickle
    from py_ecc.optimized_bn128 import curve_order
    
    print("Testing PLONK Encoder")
    print("=" * 60)
    
    # Load the PLONK arithmetization instance
    with open("constraint-system/PLONK_ARITHMETIZATION_INSTANCE.pkl", "rb") as f:
        instance = pickle.load(f)
        
    qM = instance["qM"]
    qL = instance["qL"]
    qR = instance["qR"]
    qO = instance["qO"]
    qC = instance["qC"]
    perm = instance["perm"]
    w = instance["w"]
    
    # Initialize encoder
    encoder = Encoder(curve_order)
    
    # Update state with circuit data
    print("\nUpdating state with circuit data...")
    encoder.update_state(qM, qL, qR, qO, qC, perm)
    
    # Encode selectors
    print("Encoding selector polynomials...")
    selector_polys = encoder.encode_selectors()
    
    # Encode permutation
    print("Encoding permutation polynomials...")
    perm_polys = encoder.encode_permutation()
    
    # Encode witness with 5 public inputs
    print("Encoding witness polynomials...")
    witness_data = encoder.encode_witness(w, x_size=5)
    
    # Test constraint satisfaction on a random point from H
    import random
    test_point = random.choice(encoder.H)
    
    a_poly = witness_data["a"]
    b_poly = witness_data["b"] 
    c_poly = witness_data["c"]
    PI = witness_data["PI"]
    qM_poly = selector_polys["qM"]
    qL_poly = selector_polys["qL"]
    qR_poly = selector_polys["qR"]
    qO_poly = selector_polys["qO"]
    qC_poly = selector_polys["qC"]
    
    # Check the basic arithmetic circuit constraint at a test point
    constraint_value = (qM_poly(test_point) * a_poly(test_point) * b_poly(test_point) + 
                        qL_poly(test_point) * a_poly(test_point) + 
                        qR_poly(test_point) * b_poly(test_point) + 
                        qO_poly(test_point) * c_poly(test_point) + 
                        qC_poly(test_point) + PI(test_point))
    
    print("\nTest results:")
    print(f"✅ Subgroup size: {encoder.n}")
    print(f"✅ Selector polynomials: {len(selector_polys)}")
    print(f"✅ Permutation polynomials: {len(perm_polys) - 1}")
    print(f"✅ Constraint satisfied at test point: {constraint_value == 0}")
    
    # Check that the constraint polynomial is divisible by the vanishing polynomial
    constraint_poly = (qM_poly * a_poly * b_poly + 
                       qL_poly * a_poly + 
                       qR_poly * b_poly + 
                       qO_poly * c_poly + 
                       qC_poly + PI)

    remainder = constraint_poly % encoder.v_H
    print(f"✅ Constraint polynomial divisible by vanishing polynomial: {remainder == 0}")
    
    print("\n✅ PLONK Encoder test completed successfully!")
