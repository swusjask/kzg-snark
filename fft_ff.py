from sage.all import *

def fft_ff(coeffs, w, F):
    """
    Compute the FFT of coeffs in a finite field using powers of w.
    
    Args:
        coeffs: List of field elements
        w: Primitive n-th root of unity in the field
        F: The finite field
    """
    n = len(coeffs)
    if n == 1:
        return coeffs
    
    # Split into even and odd indices
    even = coeffs[0::2]
    odd = coeffs[1::2]
    
    # Recursive FFT on the halves
    w_squared = w**2
    even_fft = fft_ff(even, w_squared, F)
    odd_fft = fft_ff(odd, w_squared, F)
    
    # Combine results
    result = [F(0)] * n
    w_power = F(1)
    
    for i in range(n//2):
        result[i] = even_fft[i] + w_power * odd_fft[i]
        result[i + n//2] = even_fft[i] - w_power * odd_fft[i]
        w_power *= w
    
    return result

def ifft_ff(values, w, F):
    """
    Compute the inverse FFT of values in a finite field.
    
    Args:
        values: List of field elements
        w: Primitive n-th root of unity in the field
        F: The finite field
    """
    n = len(values)
    # For inverse FFT, we use w^(-1) and divide results by n
    w_inv = w**(-1)
    result = fft_ff(values, w_inv, F)
    
    # Divide by n (multiply by n^(-1) in the F)
    n_inv = F(n)**(-1)
    return [x * n_inv for x in result]

def fft_ff_interpolation(values, g, F):
    """
    Compute a polynomial interpolation using FFT in a finite field.
    
    Args:
        values: List of field elements to interpolate (length must be power of 2)
        g: Generator element of order 2^k in the finite field
        F: The finite field
    """
    n = len(values)
    # Ensure n is a power of 2
    assert (n & (n-1)) == 0, "Length of values must be a power of 2"
    
    # Check that g has sufficient order
    order = g.multiplicative_order()
    assert order >= n, f"Order of g ({order}) must be at least n ({n})"
    
    # Compute coefficients using inverse FFT
    coeffs = ifft_ff(values, g, F)
    
    # Create the polynomial
    R = PolynomialRing(F, 'X')
    return R(coeffs)