from sage.all import GF, PolynomialRing


class KZG:
    """
    Kate-Zaverucha-Goldberg polynomial commitment scheme implementation.

    This class implements the KZG polynomial commitment scheme, which allows:
    1. Committing to polynomials with a single group element
    2. Creating succinct proofs for polynomial evaluations
    3. Verifying evaluations with efficient pairing checks
    4. Batching multiple polynomial evaluations for efficiency
    
    The implementation follows the formal definition in bab3.tex, supporting
    multi-polynomial commitments and batch verification optimizations.
    """

    def __init__(self, curve_type="bn254"):
        """
        Initialize the KZG polynomial commitment scheme with a pairing-friendly curve.

        Args:
            curve_type: Type of curve to use ('bn254' or 'bls12_381')
        """
        # Set up the curve operations based on the specified curve type
        if curve_type == "bn254":
            from py_ecc.optimized_bn128 import (
                G1, G2, multiply, add, curve_order, pairing, 
                neg, Z1, Z2, eq
            )
        elif curve_type == "bls12_381":
            from py_ecc.optimized_bls12_381 import (
                G1, G2, multiply, add, curve_order, pairing, 
                neg, Z1, Z2, eq
            )
        else:
            raise ValueError(f"Unsupported curve type: {curve_type}")

        # Store curve operations
        self.G1 = G1
        self.G2 = G2
        self.Z1 = Z1  # Zero point in G1
        self.Z2 = Z2  # Zero point in G2
        self.multiply = multiply
        self.add = add
        self.neg = neg
        self.pairing = pairing
        self.eq = eq
        self.curve_order = curve_order

        # Set up the finite field and polynomial ring
        self.Fq = GF(curve_order)
        self.R = PolynomialRing(self.Fq, "X")
        self.X = self.R.gen()

    def setup(self, max_degree):
        """
        Generate KZG structured reference string (commitment and verification keys).

        Args:
            max_degree: Maximum polynomial degree supported by this setup

        Returns:
            tuple: (commitment_key, verification_key)
        """
        # Sample a random secret tau ∈ Fq
        tau = self.Fq.random_element()

        # Generate commitment key: [G₁, τG₁, τ²G₁, ..., τᵈG₁]
        powers_of_tau_G1 = [self.G1]
        for i in range(1, max_degree + 1):
            powers_of_tau_G1.append(self.multiply(self.G1, int(tau**i)))

        # Generate verification key: τG₂
        tau_G2 = self.multiply(self.G2, int(tau))

        # Return the key pair
        return (powers_of_tau_G1, tau_G2)

    def commit(self, ck, polynomials):
        """
        Commit to a list of polynomials using the KZG scheme.

        Args:
            ck: Commitment key (powers of tau in G1)
            polynomials: List of polynomials to commit to

        Returns:
            List of KZG commitments (G1 points)
        """
        # Ensure all inputs are SageMath polynomials
        sage_polynomials = []
        for poly in polynomials:
            if isinstance(poly, list):
                sage_polynomials.append(self.R(poly))
            else:
                sage_polynomials.append(poly)

        max_degree = len(ck) - 1
        commitments = []

        for poly in sage_polynomials:
            if poly.degree() > max_degree:
                raise ValueError(
                    f"Polynomial degree {poly.degree()} exceeds maximum allowed degree {max_degree}"
                )

            # Compute the commitment: p(τ)·G₁ = Σᵢ pᵢ·(τⁱG₁)
            commitment = self.Z1  # Zero point in G1
            coeffs = poly.list()  # Get coefficients [p₀, p₁, ...]
            
            for i, coeff in enumerate(coeffs):
                if coeff == 0:
                    continue
                term = self.multiply(ck[i], int(coeff))
                commitment = self.add(commitment, term)

            commitments.append(commitment)

        return commitments

    def open(self, ck, polynomials, z, xi):
        """
        Create an evaluation proof for polynomials at a point.

        Args:
            ck: Commitment key
            polynomials: List of polynomials to open
            z: Evaluation point
            xi: Opening challenge for batch combining

        Returns:
            Evaluation proof π
        """
        # Ensure all inputs are SageMath polynomials
        sage_polynomials = []
        for poly in polynomials:
            if isinstance(poly, list):
                sage_polynomials.append(self.R(poly))
            else:
                sage_polynomials.append(poly)

        # Convert z and challenge to field elements
        z = self.Fq(z)
        xi = self.Fq(xi)

        # Compute the batch polynomial: p(X) = Σᵢ ξⁱ·pᵢ(X)
        combined_poly = self.R(0)  # Zero polynomial
        for i, poly in enumerate(sage_polynomials):
            combined_poly += xi ** (i + 1) * poly

        # Compute witness polynomial w(X) = (p(X) - p(z))/(X - z)
        X = self.X
        witness_poly = (combined_poly - combined_poly(z)) // (X - z)

        # Commit to the witness polynomial
        proof = self.commit(ck, [witness_poly])[0]

        return proof

    def check(self, rk, commitments, z, evaluations, proof, xi):
        """
        Verify polynomial evaluations with KZG proof.

        Args:
            rk: Verification key (τG₂)
            commitments: List of polynomial commitments
            z: Evaluation point
            evaluations: List of claimed evaluations at point z
            proof: Evaluation proof π
            xi: Opening challenge used in open

        Returns:
            bool: True if verification passes, False otherwise
        """
        # Extract the verification key component
        tau_G2 = rk

        # Convert z and challenge to field elements
        z = self.Fq(z)
        xi = self.Fq(xi)

        # Compute the batch commitment: C = Σᵢ ξⁱ·Cᵢ
        combined_commitment = self.Z1
        for i, comm in enumerate(commitments):
            challenge_power = int(xi ** (i + 1))
            term = self.multiply(comm, challenge_power)
            combined_commitment = self.add(combined_commitment, term)

        # Compute the batch evaluation: v = Σᵢ ξⁱ·vᵢ
        combined_evaluation = self.Fq(0)
        for i, eval_i in enumerate(evaluations):
            combined_evaluation += xi ** (i + 1) * self.Fq(eval_i)

        # Convert combined evaluation to integer for curve operations
        eval_int = int(combined_evaluation)
        z_int = int(z)

        # Compute C - v·G₁
        v_G1 = self.multiply(self.G1, eval_int)
        C_minus_v = self.add(combined_commitment, self.neg(v_G1))

        # Compute τG₂ - z·G₂
        z_G2 = self.multiply(self.G2, z_int)
        tauG2_minus_z = self.add(tau_G2, self.neg(z_G2))

        # Check the pairing equation: e(C - v·G₁, G₂) = e(π, τG₂ - z·G₂)
        left_pairing = self.pairing(self.G2, C_minus_v)
        right_pairing = self.pairing(tauG2_minus_z, proof)

        return left_pairing == right_pairing

    def batch_check(self, rk, commitments_list, z_list, evaluations_list, proof_list, xi_list, r=None):
        """
        Batch verify multiple polynomial evaluation proofs with a single pairing equation.

        This optimization reduces verification cost from 2n pairings to just 2 pairings
        regardless of the number of verification instances.

        Args:
            rk: Verification key (must be the same for all verifications)
            commitments_list: List of lists of commitments
            z_list: List of evaluation points
            evaluations_list: List of lists of evaluations
            proof_list: List of evaluation proofs
            xi_list: List of opening challenges
            r: Optional random field element for batching (if not provided, a new one is sampled)

        Returns:
            bool: True if all verifications pass, False otherwise
        """
        # Extract the verification key component
        tau_G2 = rk

        # Sample random field element r for batching
        if r is None:
            r = self.Fq.random_element()

        # Initialize accumulators for the batched equation
        left_acc = self.Z1  # Zero point in G1
        right_acc = self.Z1  # Zero point in G1

        # Process each verification instance
        for i, (commitments, z, evaluations, proof, xi) in enumerate(
            zip(commitments_list, z_list, evaluations_list, proof_list, xi_list)
        ):
            # Convert z and challenge to field elements
            z = self.Fq(z)
            xi = self.Fq(xi)

            # Compute the batch commitment and evaluation
            combined_commitment = self.Z1
            combined_evaluation = self.Fq(0)
            
            for j, comm in enumerate(commitments):
                xi_power = xi ** (j + 1)
                combined_commitment = self.add(
                    combined_commitment, self.multiply(comm, int(xi_power))
                )
                combined_evaluation += xi_power * self.Fq(evaluations[j])

            # Convert to integers for curve operations
            eval_int = int(combined_evaluation)
            z_int = int(z)

            # Transform the verification equation
            # From: e(C - v·G₁, G₂) = e(π, τG₂ - z·G₂)
            # To:   e(C - v·G₁ + z·π, G₂) = e(π, τG₂)
            v_G1 = self.multiply(self.G1, eval_int)
            C_minus_v = self.add(combined_commitment, self.neg(v_G1))
            z_pi = self.multiply(proof, z_int)
            term_left = self.add(C_minus_v, z_pi)

            # Apply random power r^i to this verification instance
            r_power = int(r ** (i + 1))  # Use r^(i+1) for security
            term_left = self.multiply(term_left, r_power)
            term_right = self.multiply(proof, r_power)

            # Accumulate terms for the batched equation
            left_acc = self.add(left_acc, term_left)
            right_acc = self.add(right_acc, term_right)

        # Check the batched equation:
        # e(∑(r^i(C_i - v_i·G₁ + z_i·π_i)), G₂) = e(∑(r^i·π_i), τG₂)
        left_pairing = self.pairing(self.G2, left_acc)
        right_pairing = self.pairing(tau_G2, right_acc)

        return left_pairing == right_pairing


if __name__ == "__main__":
    print("Testing KZG Polynomial Commitment Scheme")
    print("=" * 60)

    # Initialize the scheme
    kzg = KZG(curve_type="bn254")

    # Setup with max degree 5
    max_degree = 5
    ck, rk = kzg.setup(max_degree)

    # Define polynomials for testing
    R = kzg.R
    X = kzg.X
    F = kzg.Fq

    # Define multiple lists of polynomials for batch testing
    poly_list = [
        [1 + 2 * X + 3 * X**2, 4 + 5 * X**3],  # list 1
        [7 - 2 * X**2 + X**3, 3 + 4 * X + 2 * X**2],  # list 2
        [2 * X + 5 * X**2, 1 + X + X**2 + X**3],  # list 3
    ]

    # Commit to each list
    commitment_list = [kzg.commit(ck, polys) for polys in poly_list]

    # Generate random evaluation points
    z_values = [F.random_element() for _ in range(len(poly_list))]

    # Generate random opening challenges
    xi_values = [F.random_element() for _ in range(len(poly_list))]

    # Compute evaluations for verification
    evaluations_list = [
        [int(poly(z)) for poly in polys] 
        for polys, z in zip(poly_list, z_values)
    ]
    
    # Generate proofs for each list
    proof_list = [
        kzg.open(ck, polys, z, xi) 
        for polys, z, xi in zip(poly_list, z_values, xi_values)
    ]

    print("\nVerifying each commitment list individually:")
    individual_results = []
    for i, (commitments, z, evaluations, proof, xi) in enumerate(
        zip(commitment_list, z_values, evaluations_list, proof_list, xi_values)
    ):
        result = kzg.check(rk, commitments, z, evaluations, proof, xi)
        individual_results.append(result)
        print(f"{'✅' if result else '❌'} List {i+1}: {result}")

    print("\nVerifying all lists with batch verification:")
    batch_result = kzg.batch_check(
        rk, commitment_list, z_values, evaluations_list, proof_list, xi_values
    )
    print(f"{'✅' if batch_result else '❌'} Batch verification: {batch_result}")

    # Ensure batch verification matches individual verifications
    expected_batch_result = all(individual_results)
    match_result = batch_result == expected_batch_result
    print(
        f"{'✅' if match_result else '❌'} Batch verification matches individual results: {match_result}"
    )

    # Test with an invalid proof
    print("\nTesting with an invalid evaluation:")

    # Modify one evaluation to be incorrect
    if len(evaluations_list) > 0 and len(evaluations_list[0]) > 0:
        original_value = evaluations_list[0][0]
        evaluations_list[0][0] = (original_value + 1) % kzg.curve_order

        # Verify again individually
        individual_result = kzg.check(
            rk,
            commitment_list[0],
            z_values[0],
            evaluations_list[0],
            proof_list[0],
            xi_values[0],
        )
        print(f"{'✅' if not individual_result else '❌'} Modified list verification (should fail): {individual_result}")

        # Verify with batch verification
        batch_result = kzg.batch_check(
            rk, commitment_list, z_values, evaluations_list, proof_list, xi_values
        )
        print(f"{'✅' if not batch_result else '❌'} Batch verification with invalid proof (should fail): {batch_result}")
