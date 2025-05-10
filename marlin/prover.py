from sage.all import prod
import sys, os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from kzg import KZG
from fft_ff import fft_ff, fft_ff_interpolation
from encoder import Encoder
from transcript import Transcript

class Prover:
    """
    Marlin zkSNARK prover implementation.
    
    The prover generates zero-knowledge proofs for an R1CS instance following
    the optimized Marlin protocol described in the paper. It uses the polynomial
    commitment scheme to create a succinct proof of knowledge of a valid witness.
    """
    
    def __init__(self, curve_type="bn254"):
        """
        Initialize the Marlin prover with a KZG polynomial commitment scheme.
        
        Args:
            curve_type: Type of curve to use ('bn254' or 'bls12_381')
        """
        self.kzg = KZG(curve_type=curve_type)
        
    def prove(self, ipk, x, w, zero_knowledge_bound=2):
        """
        Generate a zero-knowledge proof for the R1CS instance.
        
        Args:
            ipk: Index proving key from the indexer
            x: Public input vector
            w: Witness vector (private inputs)
            zero_knowledge_bound: Degree bound for randomization polynomials
            
        Returns:
            dict: The proof containing commitments and evaluations
        """
        # Extract data from index proving key
        ck = ipk["ck"]
        A, B, C = ipk["A"], ipk["B"], ipk["C"]
        polynomials = ipk["polynomials"]
        H, K = ipk["subgroups"]["H"], ipk["subgroups"]["K"]
        n, m = ipk["subgroups"]["n"], ipk["subgroups"]["m"]
        g_K = ipk["subgroups"]["g_K"]
        v_H, v_K = ipk["vanishing_polys"]["v_H"], ipk["vanishing_polys"]["v_K"]
        R = self.kzg.R
        X = self.kzg.X
        Fq = self.kzg.Fq
        
        # Create an encoder with the same field
        self.encoder = Encoder(self.kzg.curve_order)
        self.encoder.update_state(A, B, C)
        
        # Create a transcript for the Fiat-Shamir transform
        transcript = Transcript("marlin-proof", Fq)
        transcript.append_message("public-inputs", x)
        
        # Phase 1: Encode witness and linear combinations
        z = list(x) + list(w)  # Full variable assignment
        x_size = len(x)
        
        # Compute vanishing polynomials for x and w
        v_H_x = prod([X - h for h in H[:x_size]])
        v_H_w = prod([X - h for h in H[x_size:]])
        
        # Encode witness
        encoded_witness = self.encoder.encode_witness(z, x_size)
        
        # Encode linear combinations
        encoded_combinations = self.encoder.encode_linear_combinations(z)
        
        # Extract polynomials
        w_poly = encoded_witness["w_poly"]
        x_poly = encoded_witness["x_poly"]
        zA_poly = encoded_combinations["zA_poly"]
        zB_poly = encoded_combinations["zB_poly"]
        zC_poly = encoded_combinations["zC_poly"]
        
        # Add randomness for zero-knowledge (bounded independence)
        b = zero_knowledge_bound
        
        # Random polynomials of degree < b
        w_random = sum(Fq.random_element() * X**i for i in range(b))
        zA_random = sum(Fq.random_element() * X**i for i in range(b))
        zB_random = sum(Fq.random_element() * X**i for i in range(b))
        zC_random = sum(Fq.random_element() * X**i for i in range(b))
        
        # Mask the polynomials with randomness
        w_masked = w_poly + w_random * v_H_w
        zA_masked = zA_poly + zA_random * v_H
        zB_masked = zB_poly + zB_random * v_H
        zC_masked = zC_poly + zC_random * v_H
        z_masked = w_masked * v_H_x + x_poly
        
        # Compute h_0 for the first check: zA·zB - zC = h_0·v_H
        h_0 = (zA_masked * zB_masked - zC_masked) // v_H
        assert h_0 * v_H == zA_masked * zB_masked - zC_masked, "h_0 polynomial is not well-defined"
        
        # Generate random polynomial s(X) such that sum(s(H)) = 0
        s_random = sum(Fq.random_element() * X**i for i in range(2*n+b-1))
        s_sum = sum(s_random(h) for h in H)
        s = s_random - s_sum/len(H)  # Adjust to ensure sum over H is zero
        
        # First round commitments
        first_round_polys = [w_masked, zA_masked, zB_masked, zC_masked, h_0, s]
        first_round_commitments = self.kzg.commit(ck, first_round_polys)
        
        # Add commitments to transcript
        transcript.append_message("round1-commitments", first_round_commitments)
        
        # Get first round challenges
        eta_A = transcript.get_challenge("eta_A")
        eta_B = transcript.get_challenge("eta_B")
        eta_C = transcript.get_challenge("eta_C")
        alpha = transcript.get_challenge("alpha")
        
        # Ensure alpha is not in H (as required by protocol)
        while alpha in H:
            alpha = transcript.get_challenge("alpha-retry")
        
        # Compute t(X) - the combined polynomial for r_M(α, X)
        t = self._compute_t_polynomial(
            polynomials, eta_A, eta_B, eta_C, alpha, v_H, K, R
        )
        
        # Compute first sumcheck polynomial
        r = lambda a, b: self.encoder.u_H(a, b)  # Helper function for u_H
        
        # Compute the sum of s(X) + r(α, X) * (Σ η_M * z_M(X)) - t(X) * z(X) over H
        poly = R(s + r(alpha, X)*(eta_A * zA_masked + eta_B * zB_masked + eta_C * zC_masked) - t * z_masked)
        
        # Divide by vanishing polynomial v_H to get h_1 and g_1
        h_1 = poly // v_H
        g_1 = poly % v_H

        assert g_1.constant_coefficient() == 0, "Sum over H is not 0"
        g_1 = g_1 // X
        assert h_1*v_H + X * g_1 == poly, "h_1 and g_1 are not well-defined"
        
        # Second round commitments
        second_round_polys = [t, g_1, h_1]
        second_round_commitments = self.kzg.commit(ck, second_round_polys)
        
        transcript.append_message("round2-commitments", second_round_commitments)
        
        # Get second round challenge
        beta_1 = transcript.get_challenge("beta_1")
        
        # Ensure beta_1 is not in H
        while beta_1 in H:
            beta_1 = transcript.get_challenge("beta_1-retry")
        
        # Calculate a(X) and b(X) polynomials for the third sumcheck
        a, b = self._compute_a_b_polynomials(
            polynomials, eta_A, eta_B, eta_C, beta_1, alpha, v_H, R
        )
        
        # Calculate sumcheck for the value of t(β₁)
        t_beta1 = t(beta_1)
        
        # Calculate f_2, g_2, and h_2 polynomials for the third sumcheck
        f_2 = self._compute_f2_polynomial(
            polynomials, eta_A, eta_B, eta_C, beta_1, alpha, v_H, v_K, K, g_K, Fq, R
        )
        
        # Verify f_2(0) = t(β₁)/m
        assert f_2.constant_coefficient() == t_beta1 / m, "f_2 polynomial is incorrect"

        # Compute g_2 and h_2
        g_2 = f_2 // X
        h_2 = (a - b*f_2) // v_K
        assert h_2 * v_K == a - b * (X * g_2 + t_beta1 / m), "h_2 and g_2 are not well-defined"
        
        # Third round commitments
        third_round_polys = [g_2, h_2]
        third_round_commitments = self.kzg.commit(ck, third_round_polys)
        
        transcript.append_message("round3-commitments", third_round_commitments)
        
        # Get third round challenge
        beta_2 = transcript.get_challenge("beta_2")
        
        # Evaluate polynomials at the challenge points
        polys_beta1 = [w_masked, zA_masked, zB_masked, zC_masked, h_0, s, t, g_1, h_1]
        evals_beta1 = [p(beta_1) for p in polys_beta1]
        
        # Polynomials for beta_2 evaluation
        indexer_poly_list = []
        for matrix in ["A", "B", "C"]:
            for poly_type in ["row", "col", "val"]:
                key = f"{poly_type}_{matrix}"
                indexer_poly_list.append(polynomials[key])
                
        polys_beta2 = [g_2, h_2] + indexer_poly_list
        evals_beta2 = [p(beta_2) for p in polys_beta2]
        
        # Add evaluations to transcript
        transcript.append_message("evaluations-beta1", evals_beta1)
        transcript.append_message("evaluations-beta2", evals_beta2)
        
        # Get opening challenges for batch verification
        xi_1 = transcript.get_challenge("xi_1")
        xi_2 = transcript.get_challenge("xi_2")
        
        # Generate KZG batch proofs
        proof_beta1 = self.kzg.open(ck, polys_beta1, beta_1, xi_1)
        proof_beta2 = self.kzg.open(ck, polys_beta2, beta_2, xi_2)
        
        # Assemble the final proof
        proof = {
            "commitments": {
                "first_round": first_round_commitments,
                "second_round": second_round_commitments,
                "third_round": third_round_commitments
            },
            "evaluations": {
                "beta1": evals_beta1,
                "beta2": evals_beta2
            },
            "kzg_proofs": {
                "beta1": proof_beta1,
                "beta2": proof_beta2
            }
        }
        
        return proof
        
    def _compute_t_polynomial(self, polynomials, eta_A, eta_B, eta_C, alpha, v_H, K, R):
        """
        Compute the t(X) polynomial used in the Marlin protocol.
        
        This implements the optimized version from the Marlin paper.
        
        Args:
            polynomials: Dictionary of indexed polynomials
            eta_A, eta_B, eta_C: Challenge coefficients for matrices
            alpha: Challenge point
            v_H: Vanishing polynomial for domain H
            K: Domain K
            R: Polynomial ring
            
        Returns:
            t_poly: The combined polynomial t(X)
        """
        # Extract row, col, val polynomials for each matrix
        row_A = polynomials["row_A"]
        col_A = polynomials["col_A"]
        val_A = polynomials["val_A"]
        
        row_B = polynomials["row_B"]
        col_B = polynomials["col_B"]
        val_B = polynomials["val_B"]
        
        row_C = polynomials["row_C"]
        col_C = polynomials["col_C"] 
        val_C = polynomials["val_C"]
        
        X = R.gen()
        t_poly = R(0)
        
        # Compute t(X) = Σ η_M · Σ_κ∈K [v_H(X)v_H(α)·val_M*(κ) / ((X - row_M*(κ))(α - col_M*(κ)))]
        for kappa in K:
            # Term for matrix A
            denom_A = (X - row_A(kappa)) * (alpha - col_A(kappa))
            if denom_A != 0:
                term_A = (v_H(X) * v_H(alpha) * val_A(kappa)) / denom_A
                t_poly += eta_A * term_A
            
            # Term for matrix B
            denom_B = (X - row_B(kappa)) * (alpha - col_B(kappa))
            if denom_B != 0:
                term_B = (v_H(X) * v_H(alpha) * val_B(kappa)) / denom_B
                t_poly += eta_B * term_B
            
            # Term for matrix C
            denom_C = (X - row_C(kappa)) * (alpha - col_C(kappa))
            if denom_C != 0:
                term_C = (v_H(X) * v_H(alpha) * val_C(kappa)) / denom_C
                t_poly += eta_C * term_C
        
        return R(t_poly)
        
    def _compute_a_b_polynomials(self, polynomials, eta_A, eta_B, eta_C, beta_1, alpha, v_H, R):
        """
        Compute the a(X) and b(X) polynomials for the second sumcheck in Marlin.
        
        Args:
            polynomials: Dictionary of indexed polynomials
            eta_A, eta_B, eta_C: Challenge coefficients for matrices
            beta_1, alpha: Challenge points
            v_H: Vanishing polynomial for domain H
            R: Polynomial ring
            
        Returns:
            tuple: (a_poly, b_poly)
        """
        # Extract row, col, val polynomials for each matrix
        row_A = polynomials["row_A"]
        col_A = polynomials["col_A"]
        val_A = polynomials["val_A"]
        
        row_B = polynomials["row_B"]
        col_B = polynomials["col_B"]
        val_B = polynomials["val_B"]
        
        row_C = polynomials["row_C"]
        col_C = polynomials["col_C"]
        val_C = polynomials["val_C"]
        
        a = R(0)
        b = R(1)  # Start with 1 for the product
        
        # Process each matrix
        for matrix_idx, (eta, row, col, val) in enumerate([
            (eta_A, row_A, col_A, val_A),
            (eta_B, row_B, col_B, val_B),
            (eta_C, row_C, col_C, val_C)
        ]):
            # Calculate the product term for other matrices
            other_product = R(1)
            for other_idx, (other_row, other_col) in enumerate([
                (row_A, col_A), (row_B, col_B), (row_C, col_C)
            ]):
                if other_idx != matrix_idx:
                    other_product *= (beta_1 - other_row) * (alpha - other_col)
            
            # Add term to a(X)
            a += eta * v_H(beta_1) * v_H(alpha) * val * other_product
            
            # Update b(X) with this matrix's factors
            b *= (beta_1 - row) * (alpha - col)
        
        return a, b

    def _compute_f2_polynomial(self, polynomials, eta_A, eta_B, eta_C, beta_1, alpha, v_H, v_K, K, g_K, Fq, R):
        """
        Compute the f2(X) polynomial for the second sumcheck in Marlin.
        
        Args:
            polynomials: Dictionary of indexed polynomials
            eta_A, eta_B, eta_C: Challenge coefficients for matrices
            beta_1, alpha: Challenge points
            v_H, v_K: Vanishing polynomials for domains H and K
            K: Domain K (list of field elements)
            g_K: Generator for domain K
            Fq: Finite field
            R: Polynomial ring
            
        Returns:
            f2_poly: Polynomial f₂(X)
        """
        # Extract row, col, val polynomials for each matrix
        row_A = polynomials["row_A"] 
        col_A = polynomials["col_A"]
        val_A = polynomials["val_A"]
        
        row_B = polynomials["row_B"]
        col_B = polynomials["col_B"]
        val_B = polynomials["val_B"]
        
        row_C = polynomials["row_C"]
        col_C = polynomials["col_C"]
        val_C = polynomials["val_C"]
        
        # Pre-compute v_H values
        v_H_beta1 = v_H(beta_1)
        v_H_alpha = v_H(alpha)
        
        # Pre-compute evaluations at points in K for efficiency
        row_A_evals = fft_ff(list(row_A), g_K, Fq)
        col_A_evals = fft_ff(list(col_A), g_K, Fq)
        val_A_evals = fft_ff(list(val_A), g_K, Fq)
        
        row_B_evals = fft_ff(list(row_B), g_K, Fq)
        col_B_evals = fft_ff(list(col_B), g_K, Fq)
        val_B_evals = fft_ff(list(val_B), g_K, Fq)
        
        row_C_evals = fft_ff(list(row_C), g_K, Fq)
        col_C_evals = fft_ff(list(col_C), g_K, Fq)
        val_C_evals = fft_ff(list(val_C), g_K, Fq)
        
        # Compute f2(κ) for each κ ∈ K
        f2_evals = []
        
        for i in range(len(K)):
            # Calculate denominators for each matrix
            denom_A = (beta_1 - row_A_evals[i]) * (alpha - col_A_evals[i])
            denom_B = (beta_1 - row_B_evals[i]) * (alpha - col_B_evals[i])
            denom_C = (beta_1 - row_C_evals[i]) * (alpha - col_C_evals[i])
            
            # Calculate individual terms with safeguards for division by zero
            term_A = v_H_beta1 * v_H_alpha * val_A_evals[i] / denom_A if denom_A != 0 else 0
            term_B = v_H_beta1 * v_H_alpha * val_B_evals[i] / denom_B if denom_B != 0 else 0
            term_C = v_H_beta1 * v_H_alpha * val_C_evals[i] / denom_C if denom_C != 0 else 0
            
            # Combine terms with challenge values
            f2_evals.append(eta_A * term_A + eta_B * term_B + eta_C * term_C)
        
        # Interpolate to get the polynomial
        f_2 = fft_ff_interpolation(f2_evals, g_K, Fq)

        return f_2
    
if __name__ == "__main__":
    import pickle
    from indexer import Indexer
    
    print("Testing Marlin Prover")
    print("=" * 60)
    
    # Load test instance
    with open("R1CS_INSTANCE.pkl", "rb") as f:
        RICS_INSTANCE = pickle.load(f)
    A, B, C, z = RICS_INSTANCE["A"], RICS_INSTANCE["B"], RICS_INSTANCE["C"], RICS_INSTANCE["z"]
    
    # Define public input and witness
    x_size = 5  # Adjust based on your test instance
    x = z[:x_size]
    w = z[x_size:]
    
    # Initialize the indexer and prover
    indexer = Indexer(curve_type="bn254")
    prover = Prover(curve_type="bn254")
    
    # Determine maximum degree needed
    max_degree = 200  # Adjust based on instance complexity
    
    # Preprocess the constraint system
    print("\nPreprocessing R1CS constraint system...")
    ipk, ivk = indexer.preprocess(A, B, C, max_degree)
    
    # Generate the proof
    print("\nGenerating Marlin proof...")
    proof = prover.prove(ipk, x, w)
    
    # Print proof statistics
    print("\nProof generated successfully!")
    print("\nProof components:")
    print(f"✅ First round commitments: {len(proof['commitments']['first_round'])}")
    print(f"✅ Second round commitments: {len(proof['commitments']['second_round'])}")
    print(f"✅ Third round commitments: {len(proof['commitments']['third_round'])}")
    print(f"✅ Beta1 evaluations: {len(proof['evaluations']['beta1'])}")
    print(f"✅ Beta2 evaluations: {len(proof['evaluations']['beta2'])}")
    print(f"✅ KZG proofs: 2")
    
    print("\n✅ Marlin Prover test completed successfully!")
