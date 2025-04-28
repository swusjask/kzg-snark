import sys

sys.path.insert(0, '/mnt/d/Kuliah/ctf/research/kzg-snark')
from kzg import KZG
from encoder import Encoder

class Indexer:
    """
    PLONK zkSNARK indexer implementation.
    
    The indexer handles the preprocessing phase of the PLONK protocol. It:
    1. Takes circuit description (selector values and permutation)
    2. Encodes the circuit into polynomials using the PLONK encoder
    3. Commits to these polynomials using KZG commitments
    4. Outputs an index proving key (ipk) and index verification key (ivk)
    
    This follows the preprocessing phase described in section 8 of the PLONK paper.
    """
    
    def __init__(self, curve_type="bn254"):
        """
        Initialize the PLONK indexer with a KZG polynomial commitment scheme.
        
        Args:
            curve_type: Type of curve to use ('bn254' or 'bls12_381')
        """
        self.kzg = KZG(curve_type=curve_type)
        self.encoder = Encoder(self.kzg.curve_order)
    
    def preprocess(self, qM, qL, qR, qO, qC, perm, max_degree):
        """
        Preprocess the PLONK circuit to produce the index keys.
        
        Args:
            qM, qL, qR, qO, qC: Selector values at points in H
            perm: Permutation list of length 3n
            max_degree: Maximum degree supported by the SRS
            
        Returns:
            tuple: (ipk, ivk) - index proving key and index verification key
        """
        # Setup KZG commitment scheme
        ck, rk = self.kzg.setup(max_degree)
        
        # Update encoder state with circuit description
        self.encoder.update_state(qM, qL, qR, qO, qC, perm)
        
        # Encode circuit into polynomials
        selector_polys = self.encoder.encode_selectors()
        permutation_polys = self.encoder.encode_permutation()
        
        # Extract polynomials in order
        indexer_polys = [
            selector_polys["qM"],
            selector_polys["qL"],
            selector_polys["qR"],
            selector_polys["qO"],
            selector_polys["qC"],
            permutation_polys["S_sigma1"],
            permutation_polys["S_sigma2"],
            permutation_polys["S_sigma3"]
        ]
        
        # Commit to the indexer polynomials
        index_commitments = self.kzg.commit(ck, indexer_polys)
        
        # Create index proving key
        ipk = {
            "ck": ck,
            "polynomials": {
                "qM": selector_polys["qM"],
                "qL": selector_polys["qL"],
                "qR": selector_polys["qR"],
                "qO": selector_polys["qO"],
                "qC": selector_polys["qC"],
                "S_sigma1": permutation_polys["S_sigma1"],
                "S_sigma2": permutation_polys["S_sigma2"],
                "S_sigma3": permutation_polys["S_sigma3"]
            },
            "commitments": {
                "qM": index_commitments[0],
                "qL": index_commitments[1],
                "qR": index_commitments[2],
                "qO": index_commitments[3],
                "qC": index_commitments[4],
                "S_sigma1": index_commitments[5],
                "S_sigma2": index_commitments[6],
                "S_sigma3": index_commitments[7]
            },
            # Store additional data needed by the prover
            "subgroups": {
                "H": self.encoder.H,
                "n": self.encoder.n,
                "g": self.encoder.g,
                "k1": self.encoder.k1,
                "k2": self.encoder.k2
            },
            "vanishing_poly": self.encoder.v_H,
            "sigma_star": permutation_polys["sigma_star"],
        }
        
        # Create index verification key
        ivk = {
            "rk": rk,
            "commitments": {
                "qM": index_commitments[0],
                "qL": index_commitments[1],
                "qR": index_commitments[2],
                "qO": index_commitments[3],
                "qC": index_commitments[4],
                "S_sigma1": index_commitments[5],
                "S_sigma2": index_commitments[6],
                "S_sigma3": index_commitments[7]
            },
            "subgroups": {
                "n": self.encoder.n,
                "g": self.encoder.g,
                "k1": self.encoder.k1,
                "k2": self.encoder.k2
            }
        }
        
        return ipk, ivk


if __name__ == "__main__":
    import pickle
    
    print("Testing PLONK Indexer")
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
    
    # Initialize the indexer
    indexer = Indexer(curve_type="bn254")
    
    # Determine maximum degree needed
    # In practice, we'd calculate this from the circuit requirements
    n = len(qM)
    max_degree = n + 5  # For PLONK, max_degree depends on the highest degree polynomial in the protocol
    
    # Preprocess the circuit
    print("\nPreprocessing PLONK circuit...")
    ipk, ivk = indexer.preprocess(qM, qL, qR, qO, qC, perm, max_degree)
    
    # Print some statistics
    print(f"\nCircuit size: {n} gates")
    print(f"Number of indexer polynomials: {len(ipk['polynomials'])}")
    print(f"Number of indexer commitments: {len(ipk['commitments'])}")
    
    # Verify we can access key components
    print("\nVerifying ipk contents:")
    print(f"- Contains selector polynomials: {'qM' in ipk['polynomials']}")
    print(f"- Contains permutation polynomials: {'S_sigma1' in ipk['polynomials']}")
    print(f"- Contains KZG commitment key: {len(ipk['ck']) > 0}")
    
    print("\nVerifying ivk contents:")
    print(f"- Contains selector commitments: {'qM' in ivk['commitments']}")
    print(f"- Contains permutation commitments: {'S_sigma1' in ivk['commitments']}")
    print(f"- Contains KZG verification key: {ivk['rk'] is not None}")
    
    print("\n✅ PLONK Indexer preprocessing completed successfully!")