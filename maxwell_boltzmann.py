```python
"""
maxwell_boltzmann.py

Validates the Lambert W method for estimating the standard deviation (sigma) of a
1D Maxwell-Boltzmann velocity distribution from a single measurement (vx0, y0),
as described in the PRL paper "Single-Point Parameter Estimation for Gaussian
Distributions Using the Lambert W Function." Reproduces results from Table I for
Helium-4 at 300 K.

Usage:
    python maxwell_boltzmann.py
"""

import numpy as np
from lambert_w_sigma import compute_sigma

# Physical constants
K_B = 1.380649e-23  # Boltzmann constant (J/K)
M_HE = 6.6464764e-27  # Helium-4 mass (kg)

def mb_pdf(vx, T, m):
    """Compute 1D Maxwell-Boltzmann probability density at vx."""
    sigma = np.sqrt(K_B * T / m)
    return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-vx**2 / (2 * sigma**2))

def mb_temperature(sigma, m):
    """Compute temperature T = m * sigma^2 / k_B."""
    return (m * sigma**2) / K_B

def validate_mb():
    """Validate Lambert W method for Maxwell-Boltzmann with Helium-4 at 300 K."""
    T_true = 300.0  # Target temperature (K)
    m = M_HE
    sigma_true = np.sqrt(K_B * T_true / m)  # ~789.5 m/s

    print("Maxwell-Boltzmann Validation (Helium-4, T=300 K):")
    print("vx0 (m/s) | y_0 (s/m) | z | Calc. σ (m/s) | Calc. T (K) | Rel. Error (%)")
    print("-" * 70)

    test_cases = [
        500.0,  # vx0 = 500 m/s
        600.0,  # vx0 = 600 m/s
        700.0,  # vx0 = 700 m/s
    ]

    for vx0 in test_cases:
        y0 = mb_pdf(vx0, T_true, m)
        z = -2 * np.pi * y0**2 * vx0**2
        try:
            sigma_calc = compute_sigma(vx0, y0, mu=0.0)
            T_calc = mb_temperature(sigma_calc, m)
            rel_error = abs(T_calc - T_true) / T_true * 100
            print(f"{vx0:.0f} | {y0:.3e} | {z:.3f} | {sigma_calc:.1f} | "
                  f"{T_calc:.1f} | {rel_error:.2f}")
        except ValueError as e:
            print(f"Error for vx0={vx0:.0f} m/s: {e}")

if __name__ == "__main__":
    validate_mb()
```