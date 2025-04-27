from sage.all import GF, PolynomialRing


class KZG:
    """
    Kate-Zaverucha-Goldberg polynomial commitment scheme using SageMath's polynomial ring.

    This implementation follows the formal definition from the textbook, supporting:
    1. Committing to multiple polynomials (batch commitment)
    2. Opening multiple polynomials at a single point
    3. Succinctness and efficient verification
    4. Modular curve backend (supports both BN254 and BLS12-381)
    """

    def __init__(self, curve_type="bn254"):
        """
        Initialize the KZG polynomial commitment scheme.

        Args:
            curve_type: Type of curve to use ('bn128' or 'bls12_381')
        """
        # Set up the curve operations based on the chosen curve type
        if curve_type == "bn254":
            from py_ecc.optimized_bn128 import (
                G1,
                G2,
                multiply,
                add,
                curve_order,
                pairing,
                neg,
                Z1,
                Z2,
                eq,
            )
        elif curve_type == "bls12_381":
            from py_ecc.optimized_bls12_381 import (
                G1,
                G2,
                multiply,
                add,
                curve_order,
                pairing,
                neg,
                Z1,
                Z2,
                eq,
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
        Generate the public parameters (structured reference string).

        Args:
            max_degree: Maximum degree D supported by this setup

        Returns:
            tuple: (commitment_key, verification_key)
        """
        # Sample a random secret element x
        x = self.Fq.random_element()

        # Generate commitment key: [G₁, xG₁, x²G₁, ..., xᵈG₁]
        powers_of_x_G1 = [self.G1]
        for i in range(1, max_degree + 1):
            powers_of_x_G1.append(self.multiply(self.G1, int(x**i)))

        # Generate verification key: (bp, xG₂)
        x_G2 = self.multiply(self.G2, int(x))

        ck = powers_of_x_G1
        rk = x_G2

        return (ck, rk)

    def commit(self, ck, polynomials):
        """
        Commit to a list of polynomials.

        Args:
            ck: Commitment key
            polynomials: List of SageMath polynomials

        Returns:
            List of commitments (G1 points)
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

            # Compute the commitment: p(x)·G₁ = p₀·G₁ + p₁·(xG₁) + ... + pₙ·(xⁿG₁)
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
        Create an evaluation proof for multiple polynomials at a single point.

        Args:
            ck: Commitment key
            polynomials: List of SageMath polynomials
            z: Evaluation point
            xi: Opening challenge

        Returns:
            tuple: (list of evaluations, evaluation proof π)
        """
        # Ensure all inputs are SageMath polynomials
        sage_polynomials = []
        for poly in polynomials:
            if isinstance(poly, list):
                sage_polynomials.append(self.R(poly))
            else:
                sage_polynomials.append(poly)

        # Convert z and opening challenge to field elements
        z = self.Fq(z)
        xi = self.Fq(xi)

        # Compute the linear combination: p(X) = Σᵢ ξⁱ·pᵢ(X)
        combined_poly = self.R(0)  # Zero polynomial
        for i, poly in enumerate(sage_polynomials):
            combined_poly += xi ** (i + 1) * poly

        # Compute the witness polynomial for the combined polynomial
        X = self.X
        witness_poly = (combined_poly - combined_poly(z)) // (X - z)

        # Commit to the witness polynomial
        proof = self.commit(ck, [witness_poly])[0]

        return proof

    def check(self, rk, commitments, z, evaluations, proof, xi):
        """
        Verify that the commitments open to the claimed evaluations at point z.

        Args:
            rk: Verification key
            commitments: List of commitments
            z: Evaluation point
            evaluations: List of claimed evaluations
            proof: Evaluation proof π
            xi: Opening challenge used in open

        Returns:
            bool: True if verification passes, False otherwise
        """
        # Extract the verification key components
        x_G2 = rk

        # Convert z and opening challenge to field elements
        z = self.Fq(z)
        xi = self.Fq(xi)

        # Compute the linear combination of commitments: ∑ᵢ ξⁱ·Cᵢ
        combined_commitment = self.Z1
        for i, comm in enumerate(commitments):
            challenge_power = int(xi ** (i + 1))
            term = self.multiply(comm, challenge_power)
            combined_commitment = self.add(combined_commitment, term)

        # Compute the linear combination of evaluations: ∑ᵢ ξⁱ·vᵢ
        combined_evaluation = self.Fq(0)
        for i, eval_i in enumerate(evaluations):
            combined_evaluation += xi ** (i + 1) * self.Fq(eval_i)

        # Convert combined evaluation to integer for elliptic curve operations
        eval_int = int(combined_evaluation)
        z_int = int(z)

        # Compute C - v·G₁
        v_G1 = self.multiply(self.G1, eval_int)
        C_minus_v = self.add(combined_commitment, self.neg(v_G1))

        # Compute x·G₂ - z·G₂
        z_G2 = self.multiply(self.G2, z_int)
        xG2_minus_z = self.add(x_G2, self.neg(z_G2))

        # Check the pairing equation: e(C - v·G₁, G₂) = e(π, x·G₂ - z·G₂)
        left_pairing = self.pairing(self.G2, C_minus_v)
        right_pairing = self.pairing(xG2_minus_z, proof)

        return left_pairing == right_pairing

    def batch_check(
        self, rk, commitments_list, z_list, evaluations_list, proof_list, xi_list
    ):
        """
        Batch verify multiple polynomial evaluation proofs with a single pairing equation.

        This implements the optimization described in Section 3 of the text, reducing
        verification cost from 2n pairings to just 2 pairings regardless of the number
        of verification instances.

        Args:
            rk: Verification key (must be the same for all verifications)
            commitments_list: List of lists of commitments
            z_list: List of evaluation points
            evaluations_list: List of lists of evaluations
            proof_list: List of evaluation proofs
            xi_list: List of opening challenges

        Returns:
            bool: True if all verifications pass, False otherwise
        """
        # Extract the verification key components
        x_G2 = rk

        # Sample random field element r for batching
        # This is critical for security - if r is predictable, the batching isn't sound
        r = self.Fq.random_element()

        # Initialize accumulators for the batched equation
        left_acc = self.Z1  # Zero point in G1
        right_acc = self.Z1  # Zero point in G1

        # Process each verification instance
        for i, (commitments, z, evaluations, proof, xi) in enumerate(
            zip(commitments_list, z_list, evaluations_list, proof_list, xi_list)
        ):
            # Convert z and opening challenge to field elements
            z = self.Fq(z)
            xi = self.Fq(xi)

            # Compute the linear combination of commitments and evaluations
            combined_commitment = self.Z1
            combined_evaluation = self.Fq(0)
            for j in range(len(commitments)):
                xi_power = xi ** (j + 1)
                combined_commitment = self.add(
                    combined_commitment, self.multiply(commitments[j], int(xi_power))
                )
                combined_evaluation += xi_power * self.Fq(evaluations[j])

            # Convert to integers for elliptic curve operations
            eval_int = int(combined_evaluation)
            z_int = int(z)

            # Transform the verification equation:
            # From: e(C - v·G₁, G₂) = e(π, x·G₂ - z·G₂)
            # To:   e(C - v·G₁ + z·π, G₂) = e(π, x·G₂)
            v_G1 = self.multiply(self.G1, eval_int)
            C_minus_v = self.add(combined_commitment, self.neg(v_G1))
            z_pi = self.multiply(proof, z_int)
            term_left = self.add(C_minus_v, z_pi)

            # Apply random power r^i to this verification instance
            r_power = int(
                r ** (i + 1)
            )  # Starting from r^1 rather than r^0 for security
            term_left = self.multiply(term_left, r_power)
            term_right = self.multiply(proof, r_power)

            # Accumulate terms for the batched equation
            left_acc = self.add(left_acc, term_left)
            right_acc = self.add(right_acc, term_right)

        # Check the batched equation:
        # e(∑(r^i(C_i - v_i·G₁ + z_i·π_i)), G₂) = e(∑(r^i·π_i), x·G₂)
        left_pairing = self.pairing(self.G2, left_acc)
        right_pairing = self.pairing(x_G2, right_acc)

        return left_pairing == right_pairing


if __name__ == "__main__":
    print("Testing KZG Batch Verification")
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

    # Define multiple list of polynomials
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

    # Generate proofs for each list
    proofs_and_evals = [
        kzg.open(ck, poly_list, z, xi)
        for poly_list, z, xi in zip(poly_list, z_values, xi_values)
    ]

    # Extract evaluations and proofs
    evaluations_list = [[int(poly(z)) for poly in polys] for polys, z in zip(poly_list, z_values)]
    proof_list = [kzg.open(ck, polys, z, xi) for polys, z, xi in zip(poly_list, z_values, xi_values)]

    print("\nVerifying each list individually:")
    individual_results = []
    for i, (commitments, z, evaluations, proof, xi) in enumerate(
        zip(commitment_list, z_values, evaluations_list, proof_list, xi_values)
    ):
        result = kzg.check(rk, commitments, z, evaluations, proof, xi)
        individual_results.append(result)
        print(f"{'✅' if result else '❌'} List {i+1}: {result}")

    print("\nVerifying all list with batch verification:")
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
        evaluations_list[0][0] = (
            original_value + 1
        ) % kzg.curve_order  # Change to incorrect value

        # Verify again individually
        individual_result = kzg.check(
            rk,
            commitment_list[0],
            z_values[0],
            evaluations_list[0],
            proof_list[0],
            xi_values[0],
        )
        print(
            f"{'✅' if not individual_result else '❌'} Modified list 1 verification: {individual_result}"
        )

        # Verify with batch verification
        batch_result = kzg.batch_check(
            rk, commitment_list, z_values, evaluations_list, proof_list, xi_values
        )
        print(
            f"{'✅' if not batch_result else '❌'} Batch verification with invalid proof: {batch_result}"
        )
