import sys, os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from kzg import KZG
from encoder import Encoder

class Indexer:
    """
    PLONK zkSNARK indexer implementation.
    
    The indexer performs the preprocessing phase of the PLONK protocol to create:
    1. Index Proving Key (ipk) - Used by the prover to generate proofs
    2. Index Verification Key (ivk) - Used by the verifier to check proofs
    
    This preprocessing encodes the circuit description into polynomials and 
    creates commitments for the selector and permutation polynomials.
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
        
        This method encodes the circuit description (selectors and permutation)
        into polynomial form, computes the necessary commitments, and creates
        the keys needed for proving and verifying.
        
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
        
        # Organize polynomials into ordered dictionary
        indexer_polys = {
            "qM": selector_polys["qM"],
            "qL": selector_polys["qL"],
            "qR": selector_polys["qR"],
            "qO": selector_polys["qO"],
            "qC": selector_polys["qC"],
            "S_sigma1": permutation_polys["S_sigma1"],
            "S_sigma2": permutation_polys["S_sigma2"],
            "S_sigma3": permutation_polys["S_sigma3"]
        }
        
        # Create ordered list for commitment
        poly_list = [
            indexer_polys["qM"],
            indexer_polys["qL"],
            indexer_polys["qR"],
            indexer_polys["qO"],
            indexer_polys["qC"],
            indexer_polys["S_sigma1"],
            indexer_polys["S_sigma2"],
            indexer_polys["S_sigma3"]
        ]
        
        # Commit to the indexer polynomials
        commitments_list = self.kzg.commit(ck, poly_list)
        
        # Organize commitments in a dictionary
        indexer_commitments = {
            "qM": commitments_list[0],
            "qL": commitments_list[1],
            "qR": commitments_list[2],
            "qO": commitments_list[3],
            "qC": commitments_list[4],
            "S_sigma1": commitments_list[5],
            "S_sigma2": commitments_list[6],
            "S_sigma3": commitments_list[7]
        }
        
        # Create index proving key - everything the prover needs
        ipk = {
            "ck": ck,
            "polynomials": indexer_polys,
            "commitments": indexer_commitments,
            # Additional data needed by the prover
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
        
        # Create index verification key - only what the verifier needs
        ivk = {
            "rk": rk,
            "commitments": indexer_commitments,
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
    n = len(qM)
    max_degree = n + 5  # For PLONK, max_degree depends on highest degree polynomial
    
    # Preprocess the circuit
    print("\nPreprocessing PLONK circuit...")
    ipk, ivk = indexer.preprocess(qM, qL, qR, qO, qC, perm, max_degree)
    
    # Print statistics
    print("\nPreprocessing results:")
    print(f"✅ Circuit size: {n} gates")
    print(f"✅ Number of indexer polynomials: {len(ipk['polynomials'])}")
    print(f"✅ Number of indexer commitments: {len(ipk['commitments'])}")
    
    # Verify indexer output integrity
    print("\nVerifying indexer output integrity:")
    keys_match = set(ipk['polynomials'].keys()) == set(ipk['commitments'].keys())
    print(f"✅ Polynomial and commitment keys match: {keys_match}")
    
    # Verify if all required components are present
    required_components = [
        "ck" in ipk, "rk" in ivk, 
        "polynomials" in ipk, "commitments" in ivk,
        "subgroups" in ipk, "subgroups" in ivk
    ]
    all_components_present = all(required_components)
    print(f"✅ All required components present: {all_components_present}")
    
    print("\n✅ PLONK Indexer preprocessing completed successfully!")
