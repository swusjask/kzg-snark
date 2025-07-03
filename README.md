# KZG-SNARK: Zero-Knowledge Proof Implementation

Implementation of KZG polynomial commitments and zkSNARK protocols (Marlin and PLONK) in Python/SageMath.

## Features

- **KZG Polynomial Commitment Scheme**: Kate-Zaverucha-Goldberg commitments with batching
- **Marlin Protocol**: zkSNARK implementation following Appendix E of the Marlin paper
- **PLONK Protocol**: zkSNARK implementation following Section 8 of the PLONK paper
- **Test Circuits**: Pre-generated R1CS and PLONK constraint systems

## Prerequisites

- **SageMath**: Mathematical software system ([installation guide](https://www.sagemath.org/download.html))
- **py_ecc**: Elliptic curve library (`pip install py_ecc`)

## Quick Start

```bash
# Clone repository
git clone https://github.com/swusjask/kzg-snark.git
cd kzg-snark

# Install py_ecc
pip install py_ecc

# Run demonstrations
python main.py
```

The `main.py` script demonstrates:
1. KZG polynomial commitments
2. Marlin SNARK proof generation and verification
3. PLONK SNARK proof generation and verification

## Project Structure

```
kzg-snark/
├── main.py                    # Demo script
├── kzg.py                     # KZG commitments
├── fft_ff.py                  # FFT for finite fields
├── transcript.py              # Fiat-Shamir transcript
├── constraint-system/         # Test circuits
│   ├── R1CS_INSTANCE.pkl      
│   └── PLONK_ARITHMETIZATION_INSTANCE.pkl
├── marlin/                    
│   ├── encoder.py             
│   ├── indexer.py             
│   ├── prover.py              
│   └── verifier.py            
└── plonk/                     
    ├── encoder.py             
    ├── indexer.py             
    ├── prover.py              
    └── verifier.py            
```

## License

MIT License - see LICENSE file for details

## References

- [KZG10] Kate, Zaverucha, Goldberg. "Constant-Size Commitments to Polynomials and Their Applications." ASIACRYPT 2010
- [CHMMVW20] Chiesa et al. "Marlin: Preprocessing zkSNARKs with Universal and Updatable SRS." EUROCRYPT 2020
- [GWC19] Gabizon, Williamson, Ciobotaru. "PLONK: Permutations over Lagrange-bases for Oecumenical Noninteractive arguments of Knowledge." ePrint 2019
