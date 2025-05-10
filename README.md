# KZG-SNARK: Zero-Knowledge Proof Implementation

This repository contains implementations of KZG polynomial commitments and two popular SNARK (Succinct Non-interactive ARguments of Knowledge) protocols: Marlin and PLONK. These implementations are designed to demonstrate the core cryptographic components used in modern zero-knowledge proof systems.

## Features

- **KZG Polynomial Commitment Scheme**: Implementation of Kate-Zaverucha-Goldberg polynomial commitments with batching optimizations
- **Marlin Protocol**: Implementation of the Marlin zkSNARK protocol (see Appendix E of the Marlin paper)
- **PLONK Protocol**: Implementation of the PLONK zkSNARK protocol (see Section 8 of the PLONK paper)
- **Test Constraint Systems**: Pre-generated R1CS and PLONK circuit instances for testing

## Installation

### Prerequisites

The implementation requires SageMath and the py_ecc library. Follow these steps to set up your environment:

### Installing SageMath

SageMath is an open-source mathematics software system that provides the mathematical foundations needed for this project.

**On Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install sagemath
```

**On macOS (using Homebrew):**
```bash
brew install sage
```

**On Windows:**
Download and install from [SageMath's official website](https://www.sagemath.org/download.html).

### Installing py_ecc

py_ecc is a Python library for elliptic curve cryptography, which is essential for the pairing-based cryptography used in KZG commitments.

```bash
pip install py_ecc
```

### Cloning the Repository

```bash
git clone https://github.com/yourusername/kzg-snark.git
cd kzg-snark
```

## Project Structure

```
kzg-snark/
├── __init__.py
├── fft_ff.py                  # Fast Fourier Transform for finite fields
├── kzg.py                     # KZG polynomial commitment implementation
├── transcript.py              # Fiat-Shamir transcript for non-interactive proofs
├── constraint-system/         # Example constraint systems for testing
│   ├── R1CS_INSTANCE.pkl      # R1CS instance for Marlin
│   └── PLONK_ARITHMETIZATION_INSTANCE.pkl  # PLONK circuit instance
├── marlin/                    # Marlin protocol implementation
│   ├── __init__.py
│   ├── encoder.py             # Encodes R1CS into polynomial form
│   ├── indexer.py             # Preprocessing phase
│   ├── prover.py              # Proof generation
│   └── verifier.py            # Proof verification
└── plonk/                     # PLONK protocol implementation
    ├── __init__.py
    ├── encoder.py             # Encodes circuit into polynomial form
    ├── indexer.py             # Preprocessing phase
    ├── prover.py              # Proof generation
    └── verifier.py            # Proof verification
```

## Usage

### KZG Polynomial Commitments

```python
from kzg import KZG

# Initialize the commitment scheme with BN254 curve
kzg = KZG(curve_type="bn254")

# Setup with maximum degree 10
max_degree = 10
ck, rk = kzg.setup(max_degree)

# Define polynomials (coefficients in ascending order)
R = kzg.R  # Polynomial ring
X = kzg.X  # Indeterminate
polys = [1 + 2*X + 3*X**2, 4 + 5*X**3]

# Commit to polynomials
commitments = kzg.commit(ck, polys)

# Open commitments at point z with challenge xi
z = 7
xi = 42
proof = kzg.open(ck, polys, z, xi)

# Evaluate polynomials at point z
evals = [poly(z) for poly in polys]

# Verify the proof
result = kzg.check(rk, commitments, z, evals, proof, xi)
print(f"Verification result: {result}")
```

### Marlin Protocol

```python
from marlin.indexer import Indexer
from marlin.prover import Prover
from marlin.verifier import Verifier
import pickle

# Load R1CS instance
with open("constraint-system/R1CS_INSTANCE.pkl", "rb") as f:
    instance = pickle.load(f)

A, B, C, z = instance["A"], instance["B"], instance["C"], instance["z"]
x_size = 5  # Number of public inputs
x = z[:x_size]
w = z[x_size:]

# Initialize components
indexer = Indexer(curve_type="bn254")
prover = Prover(curve_type="bn254")
verifier = Verifier(curve_type="bn254")

# Preprocess the constraint system
max_degree = 200
ipk, ivk = indexer.preprocess(A, B, C, max_degree)

# Generate proof
proof = prover.prove(ipk, x, w)

# Verify proof
result = verifier.verify(ivk, x, proof)
print(f"Verification result: {result}")
```

### PLONK Protocol

```python
from plonk.indexer import Indexer
from plonk.prover import Prover
from plonk.verifier import Verifier
import pickle

# Load PLONK circuit instance
with open("constraint-system/PLONK_ARITHMETIZATION_INSTANCE.pkl", "rb") as f:
    instance = pickle.load(f)

qM = instance["qM"]
qL = instance["qL"]
qR = instance["qR"]
qO = instance["qO"]
qC = instance["qC"]
perm = instance["perm"]
w = instance["w"]

# Define public inputs and witness
x_size = 5
x = w[:x_size]
witness = w[x_size:]

# Initialize components
indexer = Indexer(curve_type="bn254")
prover = Prover(curve_type="bn254")
verifier = Verifier(curve_type="bn254")

# Preprocess the circuit
n = len(qM)
max_degree = n + 5
ipk, ivk = indexer.preprocess(qM, qL, qR, qO, qC, perm, max_degree)

# Generate proof
proof = prover.prove(ipk, x, witness)

# Verify proof
result = verifier.verify(ivk, x, proof)
print(f"Verification result: {result}")
```

## Constraint Systems

Currently, the repository includes two pre-generated constraint systems for testing:

- `R1CS_INSTANCE.pkl`: An R1CS instance compatible with the Marlin protocol
- `PLONK_ARITHMETIZATION_INSTANCE.pkl`: A PLONK circuit instance

Note that there are no provided tools to generate new constraint systems. Users who want to test with their own circuits will need to manually create the constraint system format compatible with these implementations.

## Future Work

Potential improvements for this project:

- Add tools for generating constraint systems from high-level languages

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## References

- [KZG10] Aniket Kate, Gregory M. Zaverucha, and Ian Goldberg. "Constant-Size Commitments to Polynomials and Their Applications." ASIACRYPT 2010.
- [CHMMVW20] Alessandro Chiesa, Yuncong Hu, Mary Maller, Pratyush Mishra, Noah Vesely, and Nicholas Ward. "Marlin: Preprocessing zkSNARKs with Universal and Updatable SRS." EUROCRYPT 2020.
- [GWC19] Ariel Gabizon, Zachary J. Williamson, and Oana Ciobotaru. "PLONK: Permutations over Lagrange-bases for Oecumenical Noninteractive arguments of Knowledge." ePrint 2019.
