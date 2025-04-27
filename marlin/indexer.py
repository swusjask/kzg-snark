import sys

sys.path.insert(0, '/mnt/d/Kuliah/ctf/research/kzg-snark')
from kzg import KZG
from encoder import Encoder

class Indexer:
    """
    Marlin zkSNARK indexer implementation.
    
    The indexer handles the preprocessing phase of the zkSNARK. It:
    1. Encodes the constraint system (A, B, C matrices) into polynomials
    2. Commits to these polynomials using KZG commitments
    3. Outputs an index proving key (ipk) and index verification key (ivk)
    
    This follows the "offline phase" described in section 5.3.1 and 8.1 of the Marlin paper.
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
        Preprocess the constraint system to produce the index keys.
        
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
        A_star, B_star, C_star = A.T, B.T, C.T
        for i in range(A.ncols()):
            A_star[:, i] *= self.encoder.u_H(self.encoder.H[i], self.encoder.H[i])
            B_star[:, i] *= self.encoder.u_H(self.encoder.H[i], self.encoder.H[i])
            C_star[:, i] *= self.encoder.u_H(self.encoder.H[i], self.encoder.H[i])

        self.encoder.update_state(A_star, B_star, C_star)            
        
        # Encode matrices into polynomials
        encoded_matrices = self.encoder.encode_matrices()
        
        # Extract the nine indexer polynomials
        indexer_polys = []
        for matrix in ["A", "B", "C"]:
            for poly_type in ["row", "col", "val"]:
                key = f"{poly_type}_{matrix}"
                indexer_polys.append(encoded_matrices[key])
        
        # Commit to the indexer polynomials
        index_commitments = self.kzg.commit(ck, indexer_polys)
        
        # Create index proving key
        ipk = {
            "ck": ck,
            "A": A,
            "B": B,
            "C": C,
            "polynomials": indexer_polys,
            "commitments": index_commitments,
            # Store additional data needed by the prover
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
        
        # Create index verification key
        ivk = {
            "rk": rk,
            "commitments": index_commitments,
            "subgroups": {
                "n": self.encoder.n,
                "m": self.encoder.m,
                "g_H": self.encoder.g_H,
            }
        }
        
        return ipk, ivk
    

if __name__ == "__main__":
    import pickle
    
    print("Testing Marlin Indexer")
    print("=" * 60)
    

    with open("R1CS_INSTANCE.pkl", "rb") as f:
        RICS_INSTANCE = pickle.load(f)
    A, B, C, z = RICS_INSTANCE["A"], RICS_INSTANCE["B"], RICS_INSTANCE["C"], RICS_INSTANCE["z"]
    
    # Initialize the indexer
    indexer = Indexer(curve_type="bn254")
    
    # Determine maximum degree needed
    # In practice, we'd calculate this from the matrices and protocol requirements
    max_degree = 64
    
    # Preprocess the constraint system
    print("\nPreprocessing constraint system...")
    ipk, ivk = indexer.preprocess(A, B, C, max_degree)
    
    # Print some statistics
    print(f"\nMatrix dimensions: {A.nrows()} × {A.ncols()}")
    print(f"Non-zero elements: A: {len(A.nonzero_positions())}, B: {len(B.nonzero_positions())}, C: {len(C.nonzero_positions())}")
    print(f"Subgroup sizes: H: {ipk['subgroups']['n']}, K: {ipk['subgroups']['m']}")
    print(f"Number of indexer polynomials: {len(ipk['polynomials'])}")
    print(f"Number of indexer commitments: {len(ipk['commitments'])}")
    
    print("\n✅ Indexer preprocessing completed successfully!")