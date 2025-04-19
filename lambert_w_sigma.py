```python
"""
lambert_w_sigma.py

Computes the standard deviation (sigma) of a normalized Gaussian PDF from a single
measurement (x0, y0) with known mean (mu) using the Lambert W function, as described
in the PRL paper "Single-Point Parameter Estimation for Gaussian Distributions Using
the Lambert W Function."

Usage:
    sigma = compute_sigma(x0, y0, mu)
    # Example: Trapped ion, x0=5e-9 m, y0=1.59e7 1/m, mu=0
    # Example: NMR CH3 peak, x0=chemical shift, y0=normalized intensity
"""

import numpy as np
from scipy.special import lambertw

def compute_sigma(x0, y0, mu=0.0):
    """
    Compute sigma for a normalized Gaussian PDF using the Lambert W function.

    Parameters:
    -----------
    x0 : float
        Measurement position (e.g., chemical shift in NMR, position in QHO).
    y0 : float
        Measured PDF value at x0 (e.g., normalized intensity, probability density).
    mu : float, optional
        Known mean of the Gaussian (default: 0.0).

    Returns:
    --------
    sigma : float
        Estimated standard deviation of the Gaussian.

    Raises:
    -------
    ValueError
        If inputs violate domain constraints or are invalid.
    """
    # Validate inputs
    if y0 <= 0:
        raise ValueError("y0 must be positive")
    if x0 == mu:
        raise ValueError("x0 must not equal mu")

    # Compute z = -2 * pi * y0^2 * (x0 - mu)^2
    delta = x0 - mu
    z = -2 * np.pi * y0**2 * delta**2

    # Check domain: z >= -1/e
    if z < -1 / np.e:
        raise ValueError(
            f"Input violates domain constraint: y0 * |x0 - mu| <= {1 / np.sqrt(2 * np.pi * np.e):.5f}"
        )

    # Compute sigma using Lambert W (principal branch, k=0)
    w0 = lambertw(z, k=0).real
    sigma = (1 / (y0 * np.sqrt(2 * np.pi))) * np.exp(w0 / 2)

    return sigma

if __name__ == "__main__":
    # Example: Trapped ion (QHO, 40Ca+ ion)
    x0 = 5e-9  # 5 nm
    y0 = 1.59e7  # 1.59e7 1/m
    mu = 0.0
    try:
        sigma = compute_sigma(x0, y0, mu)
        print(f"Trapped ion: sigma = {sigma*1e9:.1f} nm")
    except ValueError as e:
        print(f"Error: {e}")

    # Example: Placeholder for NMR (requires normalized data)
    # x0 = chemical_shift  # e.g., 1.2 ppm for CH3 peak
    # y0 = normalized_intensity  # e.g., from SDBS No. 1300, normalized
    # sigma = compute_sigma(x0, y0, mu=known_peak_center)

