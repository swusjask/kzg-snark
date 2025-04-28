import sys
from sage.all import prod

sys.path.insert(0, '/mnt/d/Kuliah/ctf/research/kzg-snark')
from kzg import KZG
from fft_ff import fft_ff_interpolation
from encoder import Encoder
from transcript import Transcript

class Prover:
    """
    PLONK zkSNARK prover implementation.
    
    The prover generates a zero-knowledge proof for a PLONK circuit following
    the protocol described in section 8.3 of the PLONK paper.
    """
    
    def __init__(self, curve_type="bn254"):
        """
        Initialize the PLONK prover with a KZG polynomial commitment scheme.
        
        Args:
            curve_type: Type of curve to use ('bn254' or 'bls12_381')
        """
        self.kzg = KZG(curve_type=curve_type)
        
    def prove(self, ipk, x, w):
        """
        Generate a zero-knowledge proof for the PLONK circuit.
        
        Args:
            ipk: Index proving key from the indexer
            x: Public input list
            w: Witness list (private inputs)
            
        Returns:
            dict: The proof
        """
        # Extract data from index proving key
        ck = ipk["ck"]
        polynomials = ipk["polynomials"]
        H = ipk["subgroups"]["H"]
        n = ipk["subgroups"]["n"]
        g = ipk["subgroups"]["g"]
        k1 = ipk["subgroups"]["k1"]
        k2 = ipk["subgroups"]["k2"]
        v_H = ipk["vanishing_poly"]
        sigma_star = ipk["sigma_star"]
        Fq = self.kzg.Fq
        R = self.kzg.R
        X = self.kzg.X
        
        # Create an encoder with the same field
        self.encoder = Encoder(self.kzg.curve_order)
        
        # Create a transcript for the Fiat-Shamir transform
        transcript = Transcript("plonk-proof", Fq)
        
        # Add public inputs to the transcript
        transcript.append_message("public-inputs", x)
        
        # Compute the full witness
        full_witness = x + w
        
        # Initialize the encoder with a temporary empty permutation (not needed for witness encoding)
        # We're just using the encoder to help with computing the public input polynomial
        empty_perm = [0] * (3 * n)
        empty_selectors = [Fq(0)] * n
        self.encoder.update_state(empty_selectors, empty_selectors, empty_selectors, empty_selectors, empty_selectors, empty_perm)
        
        # Compute public input polynomial
        PI = self.encoder.compute_public_input_poly(x)
        
        # ----- Round 1: Wire polynomials -----
        # Generate random blinding scalars
        b1, b2 = Fq.random_element(), Fq.random_element()
        b3, b4 = Fq.random_element(), Fq.random_element()
        b5, b6 = Fq.random_element(), Fq.random_element()
        b7, b8, b9 = Fq.random_element(), Fq.random_element(), Fq.random_element()
        
        # Split witness into a, b, c values
        a_values = full_witness[:n]
        b_values = full_witness[n:2*n]
        c_values = full_witness[2*n:3*n]
        
        # Compute wire polynomials a(X), b(X), c(X) with blinding factors
        a_poly = (b1 * X + b2) * v_H + fft_ff_interpolation(a_values, g, Fq)
        b_poly = (b3 * X + b4) * v_H + fft_ff_interpolation(b_values, g, Fq)
        c_poly = (b5 * X + b6) * v_H + fft_ff_interpolation(c_values, g, Fq)
        
        # Commit to wire polynomials
        first_round_polys = [a_poly, b_poly, c_poly]
        first_round_commitments = self.kzg.commit(ck, first_round_polys)
        a_commit, b_commit, c_commit = first_round_commitments
        
        # Add commitments to transcript
        transcript.append_message("round1-commitments", first_round_commitments)
        
        # ----- Round 2: Permutation polynomial -----
        # Get permutation challenges
        beta = transcript.get_challenge("beta")
        gamma = transcript.get_challenge("gamma")
        
        # Compute permutation polynomial z(X)
        z_poly = self._compute_permutation_polynomial(
            a_values, b_values, c_values, 
            sigma_star,
            beta, gamma, g, k1, k2, n, H, v_H, X, Fq,
            b7, b8, b9
        )
        L1 = R((X**n - 1) / (n * (X - 1)))
        assert L1 * (z_poly - 1) % v_H == 0, "z_poly does not satisfy L1 condition"

        # Commit to permutation polynomial
        z_commit = self.kzg.commit(ck, [z_poly])[0]
        
        # Add commitment to transcript
        transcript.append_message("round2-commitment", z_commit)
        
        # ----- Round 3: Quotient polynomial -----
        # Get quotient challenge
        alpha = transcript.get_challenge("alpha")
        
        # Compute quotient polynomial t(X)
        t_poly = self._compute_quotient_polynomial(
            a_poly, b_poly, c_poly, z_poly,
            polynomials["qM"], polynomials["qL"], polynomials["qR"], 
            polynomials["qO"], polynomials["qC"],
            polynomials["S_sigma1"], polynomials["S_sigma2"], polynomials["S_sigma3"],
            alpha, beta, gamma, PI, v_H, H, n, g, k1, k2, R, X
        )
        
        # Split t(X) into degree < n polynomials
        t_lo, t_mid, t_hi = self._split_quotient_polynomial(t_poly, n, X, R, Fq)
        # try testing
        # Commit to the parts of t(X)
        third_round_polys = [t_lo, t_mid, t_hi]
        third_round_commitments = self.kzg.commit(ck, third_round_polys)
        t_lo_commit, t_mid_commit, t_hi_commit = third_round_commitments
        
        # Add commitments to transcript
        transcript.append_message("round3-commitments", third_round_commitments)
        
        # ----- Round 4: Evaluation point -----
        # Get evaluation challenge
        zeta = transcript.get_challenge("zeta")
        
        # Ensure zeta is not in H
        while zeta in H:
            zeta = transcript.get_challenge("zeta-retry")
        
        # Compute opening evaluations
        a_zeta = a_poly(zeta)
        b_zeta = b_poly(zeta)
        c_zeta = c_poly(zeta)
        s_sigma1_zeta = polynomials["S_sigma1"](zeta)
        s_sigma2_zeta = polynomials["S_sigma2"](zeta)
        z_omega_zeta = z_poly(zeta * g)  # z(ζω)
        
        # Add evaluations to transcript
        transcript.append_message("round4-evaluations", [a_zeta, b_zeta, c_zeta, s_sigma1_zeta, s_sigma2_zeta, z_omega_zeta])
        
        # ----- Round 5: Opening proof -----
        # Get opening challenge
        v = transcript.get_challenge("v")
        
        # Compute linearization polynomial r(X)
        r_poly = self._compute_linearization_polynomial(
            a_zeta, b_zeta, c_zeta, s_sigma1_zeta, s_sigma2_zeta, z_omega_zeta,
            polynomials["qM"], polynomials["qL"], polynomials["qR"], 
            polynomials["qO"], polynomials["qC"], polynomials["S_sigma3"],
            z_poly, t_lo, t_mid, t_hi, alpha, beta, gamma, zeta, PI, n, k1, k2, R
        )

        assert r_poly(zeta) == 0, "r(ζ) should be zero"
        
        # First batch: Polynomials to be opened at zeta
        zeta_polys = [
            r_poly,
            a_poly, 
            b_poly, 
            c_poly,
            polynomials["S_sigma1"],
            polynomials["S_sigma2"]
        ]
        
        # Get the opening proof for polynomials at zeta
        W_z = self.kzg.open(ck, zeta_polys, zeta, v)
        W_zw = self.kzg.open(ck, [z_poly], zeta * g, v)
        
        
        # Assemble the final proof
        proof = {
            "commitments": {
                "a": a_commit,
                "b": b_commit,
                "c": c_commit,
                "z": z_commit,
                "t_lo": t_lo_commit,
                "t_mid": t_mid_commit,
                "t_hi": t_hi_commit,
                "W_z": W_z,
                "W_zw": W_zw
            },
            "evaluations": {
                "a": a_zeta,
                "b": b_zeta,
                "c": c_zeta,
                "s_sigma1": s_sigma1_zeta,
                "s_sigma2": s_sigma2_zeta,
                "z_omega": z_omega_zeta
            }
        }
        
        return proof
    
    def _compute_permutation_polynomial(self, a_values, b_values, c_values, 
                                        sigma_star,
                                        beta, gamma, g, k1, k2, n, H, v_H, X, Fq,
                                        b7, b8, b9):
        """
        Compute the permutation polynomial z(X).
        
        Following section 8.3, Round 2 of the PLONK paper.
        """
        # Add blinding factors
        z_poly = (b7 * X**2 + b8 * X + b9) * v_H
        
        # Compute the permutation accumulator at each point
        z_values = [Fq(1)]  # z(ω⁰) = 1
        
        for i in range(n - 1):
            # Compute numerator: (a(ωⁱ) + β·ωⁱ + γ)·(b(ωⁱ) + β·k₁·ωⁱ + γ)·(c(ωⁱ) + β·k₂·ωⁱ + γ)
            num = ((a_values[i] + beta * H[i] + gamma) * 
                   (b_values[i] + beta * k1 * H[i] + gamma) * 
                   (c_values[i] + beta * k2 * H[i] + gamma))
            
            # Compute denominator: (a(ωⁱ) + β·σ₁(ωⁱ) + γ)·(b(ωⁱ) + β·σ₂(ωⁱ) + γ)·(c(ωⁱ) + β·σ₃(ωⁱ) + γ)
            den = ((a_values[i] + beta * sigma_star[i] + gamma) * 
                   (b_values[i] + beta * sigma_star[i + n] + gamma) * 
                   (c_values[i] + beta * sigma_star[i + 2*n] + gamma))
            
            if den == 0:
                # This should be extremely rare and would indicate an issue with the circuit
                raise ValueError("Denominator is zero in permutation polynomial calculation")
            
            # Update the accumulator: z(ωⁱ⁺¹) = z(ωⁱ)·num/den
            z_values.append(z_values[-1] * (num / den))
        
        # Interpolate the accumulator values to get the rest of z(X)
        z_interp = fft_ff_interpolation(z_values, g, Fq)

        # Combine with the blinding part
        z_poly += z_interp
        
        return z_poly
    
    def _compute_quotient_polynomial(self, a_poly, b_poly, c_poly, z_poly,
                                     qM, qL, qR, qO, qC,
                                     S_sigma1, S_sigma2, S_sigma3,
                                     alpha, beta, gamma, PI, v_H, H, n, g, k1, k2, R, X):
        """
        Compute the quotient polynomial t(X).
        
        Following section 8.3, Round 3 of the PLONK paper.
        """
        # Term 1: Gate constraints
        term1 = R((a_poly * b_poly * qM + a_poly * qL + b_poly * qR + c_poly * qO + PI + qC) / v_H)
        
        # Term 2: Permutation constraints (part 1)
        term2 = alpha * (z_poly * (a_poly + beta * X + gamma) * 
                         (b_poly + beta * k1 * X + gamma) * 
                         (c_poly + beta * k2 * X + gamma)) / v_H
        
        # Term 3: Permutation constraints (part 2)
        z_poly_shifted = R(z_poly(g * X))  # z(ωX)
        term3 = -alpha * ((a_poly + beta * S_sigma1 + gamma) * 
                          (b_poly + beta * S_sigma2 + gamma) * 
                          (c_poly + beta * S_sigma3 + gamma) * 
                          z_poly_shifted) / v_H
        
        L1 = R((X**n - 1) / (n * (X - 1)))
        term4 = R(alpha**2 * (z_poly - 1) * L1 / v_H)
        
        # Combine all terms
        t_poly = R(term1 + term2 + term3 + term4)
        
        return t_poly
    
    def _split_quotient_polynomial(self, t_poly, n, X, R, Fq):
        """
        Split the quotient polynomial t(X) into three parts of degree < n.
        
        t(X) = t_lo(X) + X^n · t_mid(X) + X^(2n) · t_hi(X)
        
        With random blinding factors for zero-knowledge.
        """
        # Get coefficients of t_poly
        t_coeffs = t_poly.list()
        t_coeffs.extend([Fq(0)] * (3 * n - len(t_coeffs)))  # Pad to ensure 3n coefficients
        
        # Split into three parts
        t_lo_coeffs = t_coeffs[:n]
        t_mid_coeffs = t_coeffs[n:2*n]
        t_hi_coeffs = t_coeffs[2*n:]
        
        # Add random blinding factors
        b10, b11 = Fq.random_element(), Fq.random_element()
        
        # Create the polynomials with blinding
        t_lo = R(t_lo_coeffs) + b10 * X**n
        t_mid = R(t_mid_coeffs) - b10 + b11 * X**n
        t_hi = R(t_hi_coeffs) - b11

        assert t_poly == t_lo + X**n * t_mid + X**(2*n) * t_hi, "t(X) does not equal the sum of its parts"
        
        return t_lo, t_mid, t_hi
    
    def _compute_linearization_polynomial(self, a_zeta, b_zeta, c_zeta, 
                                          s_sigma1_zeta, s_sigma2_zeta, z_omega_zeta,
                                          qM, qL, qR, qO, qC, S_sigma3,
                                          z_poly, t_lo, t_mid, t_hi, alpha, beta, gamma, zeta, PI, n, k1, k2, R):
        """
        Compute the linearization polynomial r(X).
        
        Following section 8.3, Round 5 of the PLONK paper.
        """
        # Compute z_H(ζ) = ζ^n - 1
        z_H_zeta = zeta**n - 1
        
        # Compute L1(ζ)
        L1_zeta = R((z_H_zeta) / (n * (zeta - 1)))
        
        # Evaluate PI(ζ)
        PI_zeta = PI(zeta)
        
        # Compute the individual terms
        
        # Gate constraints
        term1 = a_zeta * b_zeta * qM + a_zeta * qL + b_zeta * qR + c_zeta * qO + PI_zeta + qC
        
        # Permutation constraints (part 1)
        term2 = alpha * ((a_zeta + beta * zeta + gamma) * 
                         (b_zeta + beta * k1 * zeta + gamma) * 
                         (c_zeta + beta * k2 * zeta + gamma) * z_poly)
        
        # Permutation constraints (part 2)
        term3 = -alpha * ((a_zeta + beta * s_sigma1_zeta + gamma) * 
                          (b_zeta + beta * s_sigma2_zeta + gamma) * 
                          (c_zeta + beta * S_sigma3 + gamma) * z_omega_zeta)
        
        # Copy constraints
        term4 = alpha**2 * (z_poly - 1) * L1_zeta
        
        # Compute the final linearization polynomial
        # Note we're subtracting the quotient polynomial terms
        r_poly = term1 + term2 + term3 + term4 - z_H_zeta * (
            t_lo + zeta**n * t_mid + zeta**(2*n) * t_hi
        )
        
        return r_poly
        
    def _compute_opening_proof_polynomials(self, r_poly, a_poly, b_poly, c_poly, 
                                          S_sigma1, S_sigma2, z_poly,
                                          a_zeta, b_zeta, c_zeta, 
                                          s_sigma1_zeta, s_sigma2_zeta, z_omega_zeta,
                                          zeta, v, R, X):
        """
        Compute the opening proof polynomials W_z(X) and W_zω(X).
        
        Following section A.1 of the PLONK paper (batch opening).
        """
        # Compute W_z(X) - the opening proof for polynomials evaluated at ζ
        W_z_denom = X - zeta
        W_z_num = (r_poly + 
                  v * (a_poly - a_zeta) + 
                  v**2 * (b_poly - b_zeta) + 
                  v**3 * (c_poly - c_zeta) + 
                  v**4 * (S_sigma1 - s_sigma1_zeta) + 
                  v**5 * (S_sigma2 - s_sigma2_zeta))
        
        W_z_poly = W_z_num // W_z_denom
        
        # Compute W_zω(X) - the opening proof for z(X) evaluated at ζω
        W_zw_denom = X - zeta * R.gen()**1  # ζω
        W_zw_num = z_poly - z_omega_zeta
        
        W_zw_poly = W_zw_num // W_zw_denom
        
        return W_z_poly, W_zw_poly


if __name__ == "__main__":
    import pickle
    from indexer import Indexer
    
    print("Testing PLONK Prover")
    print("=" * 60)
    
    # Load the PLONK arithmetization instance
    with open("PLONK_ARITHMETIZATION_INSTANCE.pkl", "rb") as f:
        instance = pickle.load(f)
        
    qM = instance["qM"]
    qL = instance["qL"]
    qR = instance["qR"]
    qO = instance["qO"]
    qC = instance["qC"]
    perm = instance["perm"]
    w = instance["w"]
    
    # Define public input (first few elements) and witness (remaining elements)
    x_size = 5  # Adjust based on your instance
    x = w[:x_size]
    remaining_w = w[x_size:]
    
    # Initialize the indexer and prover
    indexer = Indexer(curve_type="bn254")
    prover = Prover(curve_type="bn254")
    
    # Determine maximum degree needed
    n = len(qM)
    max_degree = n + 5  # For PLONK, max_degree depends on highest degree polynomial
    
    # Preprocess the circuit
    print("\nPreprocessing PLONK circuit...")
    ipk, ivk = indexer.preprocess(qM, qL, qR, qO, qC, perm, max_degree)
    
    # Generate the proof
    print("\nGenerating proof...")
    proof = prover.prove(ipk, x, remaining_w)
    
    # Print proof statistics
    print(f"\nProof generated successfully!")
    print(f"\nProof components:")
    print(f"- Number of commitments: {len(proof['commitments'])}")
    print(f"- Number of evaluations: {len(proof['evaluations'])}")
    print(f"- Number of challenges: {len(proof['challenges'])}")
    
    print("\n✅ PLONK Prover test completed successfully!")