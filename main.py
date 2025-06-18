#!/usr/bin/env python3
"""
KZG-SNARK Demo: KZG Commitments, Marlin, and PLONK
"""

import pickle
from kzg import KZG
from marlin.indexer import Indexer as MarlinIndexer
from marlin.prover import Prover as MarlinProver
from marlin.verifier import Verifier as MarlinVerifier
from plonk.indexer import Indexer as PlonkIndexer
from plonk.prover import Prover as PlonkProver
from plonk.verifier import Verifier as PlonkVerifier


def demo_kzg():
    print("=== KZG Polynomial Commitment Demo ===")
    
    # Initialize and setup
    kzg = KZG(curve_type="bn254")
    ck, rk = kzg.setup(max_degree=10)
    
    # Create and commit to polynomials
    R = kzg.R
    X = kzg.X
    polys = [1 + 2*X + 3*X**2, 4 + 5*X**3]
    commitments = kzg.commit(ck, polys)
    
    # Open at point z=7
    z, xi = 7, 42
    proof = kzg.open(ck, polys, z, xi)
    evals = [poly(z) for poly in polys]
    
    # Verify
    result = kzg.check(rk, commitments, z, evals, proof, xi)
    print(f"KZG verification: {'PASS' if result else 'FAIL'}\n")


def demo_marlin():
    print("=== Marlin SNARK Demo ===")
    
    # Load R1CS instance
    with open("constraint-system/R1CS_INSTANCE.pkl", "rb") as f:
        instance = pickle.load(f)
    
    A, B, C, z = instance["A"], instance["B"], instance["C"], instance["z"]
    x = z[:5]  # Public inputs
    w = z[5:]  # Private witness
    
    # Setup
    indexer = MarlinIndexer(curve_type="bn254")
    ipk, ivk = indexer.preprocess(A, B, C, max_degree=200)
    
    # Prove
    prover = MarlinProver(curve_type="bn254")
    proof = prover.prove(ipk, x, w)
    
    # Verify
    verifier = MarlinVerifier(curve_type="bn254")
    result = verifier.verify(ivk, x, proof)
    print(f"Marlin verification: {'PASS' if result else 'FAIL'}\n")


def demo_plonk():
    print("=== PLONK SNARK Demo ===")
    
    # Load PLONK circuit
    with open("constraint-system/PLONK_ARITHMETIZATION_INSTANCE.pkl", "rb") as f:
        instance = pickle.load(f)
    
    qM = instance["qM"]
    qL = instance["qL"]
    qR = instance["qR"]
    qO = instance["qO"]
    qC = instance["qC"]
    perm = instance["perm"]
    w = instance["w"]
    
    x = w[:5]  # Public inputs
    witness = w[5:]  # Private witness
    
    # Setup
    indexer = PlonkIndexer(curve_type="bn254")
    n = len(qM)
    ipk, ivk = indexer.preprocess(qM, qL, qR, qO, qC, perm, max_degree=n+5)
    
    # Prove
    prover = PlonkProver(curve_type="bn254")
    proof = prover.prove(ipk, x, witness)
    
    # Verify
    verifier = PlonkVerifier(curve_type="bn254")
    result = verifier.verify(ivk, x, proof)
    print(f"PLONK verification: {'PASS' if result else 'FAIL'}\n")


if __name__ == "__main__":
    print("Running KZG-SNARK demonstrations...\n")
    
    try:
        demo_kzg()
    except Exception as e:
        print(f"KZG demo failed: {e}\n")
    
    try:
        demo_marlin()
    except Exception as e:
        print(f"Marlin demo failed: {e}\n")
    
    try:
        demo_plonk()
    except Exception as e:
        print(f"PLONK demo failed: {e}\n")
    
    print("Demo complete!")