import sys, os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from kzg import KZG
from encoder import Encoder

class Indexer:
    """
    Marlin zkSNARK indexer implementation.
    
    The indexer performs the preprocessing phase of the Marlin protocol to create:
    1. Index Proving Key (ipk) - Used by the prover to generate proofs
    2. Index Verification Key (ivk) - Used by the verifier to check proofs
    
    This preprocessing encodes the R1CS constraint system into polynomials and 
    commits to them, generating a structured reference string specific to this instance.
    """
    
    def __init__(self, curve_type="bn254"):
        """
        Initialize the Marlin indexer with a KZG polynomial commitment scheme.
        
        Args:
            curve_type: Type of curve to use ('bn254' or 'bls12_381')
        """
        self.kzg = KZG(curve_type=curve_type)
        self.encoder = Encoder(self.kzg.curve_order)
        
    def preprocess(self, A, B, C, max_degree):
        """
        Preprocess the R1CS constraint system to produce index keys.
        
        This method encodes the R1CS matrices into polynomial form, computes
        the necessary commitments, and creates the keys needed for proving and
        verifying.
        
        Args:
            A, B, C: The R1CS constraint matrices
            max_degree: Maximum degree supported by the SRS
            
        Returns:
            tuple: (ipk, ivk) - index proving key and index verification key
        """
        # Setup KZG commitment scheme
        ck, rk = self.kzg.setup(max_degree)
        
        # Update encoder state with matrices
        self.encoder.update_state(A, B, C)
        
        # Create star matrices for more efficient R1CS representation
        A_star, B_star, C_star = A.T, B.T, C.T
        for i in range(A.ncols()):
            A_star[:, i] *= self.encoder.u_H(self.encoder.H[i], self.encoder.H[i])
            B_star[:, i] *= self.encoder.u_H(self.encoder.H[i], self.encoder.H[i])
            C_star[:, i] *= self.encoder.u_H(self.encoder.H[i], self.encoder.H[i])

        self.encoder.update_state(A_star, B_star, C_star)            
        
        # Encode matrices into polynomials
        encoded_matrices = self.encoder.encode_matrices()
        
        # Organize indexer polynomials in a dictionary
        indexer_polys = {}
        for matrix in ["A", "B", "C"]:
            for poly_type in ["row", "col", "val"]:
                key = f"{poly_type}_{matrix}"
                indexer_polys[key] = encoded_matrices[key]
        
        # Create a list version for commitment (to be consistent with paper)
        indexer_polys_list = []
        for matrix in ["A", "B", "C"]:
            for poly_type in ["row", "col", "val"]:
                key = f"{poly_type}_{matrix}"
                indexer_polys_list.append(encoded_matrices[key])
        
        # Commit to the indexer polynomials
        index_commitments = self.kzg.commit(ck, indexer_polys_list)
        
        # Organize commitments in a dictionary
        commitments = {}
        i = 0
        for matrix in ["A", "B", "C"]:
            for poly_type in ["row", "col", "val"]:
                key = f"{poly_type}_{matrix}"
                commitments[key] = index_commitments[i]
                i += 1
        
        # Create index proving key
        ipk = {
            "ck": ck,
            "A": A,
            "B": B,
            "C": C,
            "polynomials": indexer_polys,
            "commitments": commitments,
            # Additional data needed by the prover
            "subgroups": {
                "H": self.encoder.H,
                "K": self.encoder.K,
                "g_H": self.encoder.g_H,
                "g_K": self.encoder.g_K,
                "n": self.encoder.n,
                "m": self.encoder.m
            },
            "vanishing_polys": {
                "v_H": self.encoder.v_H,
                "v_K": self.encoder.v_K,
            }
        }
        
        # Create index verification key - only contains what verifier needs
        ivk = {
            "rk": rk,
            "commitments": commitments,
            "subgroups": {
                "n": self.encoder.n,
                "m": self.encoder.m,
                "g_H": self.encoder.g_H,
            },
            "vanishing_polys": {
                "v_H": self.encoder.v_H,
                "v_K": self.encoder.v_K,
            }
        }
        
        return ipk, ivk
    

if __name__ == "__main__":
    import pickle
    
    print("Testing Marlin Indexer")
    print("=" * 60)
    
    # Load test instance
    with open("R1CS_INSTANCE.pkl", "rb") as f:
        RICS_INSTANCE = pickle.load(f)
    A, B, C, z = RICS_INSTANCE["A"], RICS_INSTANCE["B"], RICS_INSTANCE["C"], RICS_INSTANCE["z"]
    
    # Initialize the indexer
    indexer = Indexer(curve_type="bn254")
    
    # Determine maximum degree for KZG setup
    max_degree = 64
    
    # Preprocess the constraint system
    print("\nPreprocessing R1CS constraint system...")
    ipk, ivk = indexer.preprocess(A, B, C, max_degree)
    
    # Print statistics
    print("\nPreprocessing results:")
    print(f"✅ Matrix dimensions: {A.nrows()} × {A.ncols()}")
    print(f"✅ Non-zero elements: A: {len(A.nonzero_positions())}, " 
          f"B: {len(B.nonzero_positions())}, C: {len(C.nonzero_positions())}")
    print(f"✅ Subgroup sizes: H: {ipk['subgroups']['n']}, K: {ipk['subgroups']['m']}")
    print(f"✅ Number of indexer polynomials: {len(ipk['polynomials'])}")
    print(f"✅ Number of indexer commitments: {len(ipk['commitments'])}")
    
    # Verify indexer output integrity
    print("\nVerifying indexer output integrity:")
    keys_match = set(ipk["polynomials"].keys()) == set(ipk["commitments"].keys())
    print(f"✅ Polynomial and commitment keys match: {keys_match}")
    
    # Verify if all required components are present
    required_components = [
        "ck" in ipk, "rk" in ivk, 
        "polynomials" in ipk, "commitments" in ivk,
        "subgroups" in ipk, "subgroups" in ivk
    ]
    all_components_present = all(required_components)
    print(f"✅ All required components present: {all_components_present}")
    
    print("\n✅ Marlin Indexer preprocessing completed successfully!")
