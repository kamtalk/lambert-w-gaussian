```python
"""
qho_validation.py

Validates the Lambert W method for estimating the standard deviation (sigma) of a
quantum harmonic oscillator (QHO) ground state from a single measurement (x0, y0),
as described in the PRL paper "Single-Point Parameter Estimation for Gaussian
Distributions Using the Lambert W Function." Reproduces results from Table I and
trapped ion (40Ca+) example.

Usage:
    python qho_validation.py
"""

import numpy as np
from lambert_w_sigma import compute_sigma

# Physical constants
HBAR = 1.0545718e-34  # Planck's constant / 2pi (J s)
M_ELECTRON = 9.1093837e-31  # Electron mass (kg)
M_CA_ION = 6.642157e-26  # 40Ca+ ion mass (kg)

def qho_pdf(x, sigma):
    """Compute QHO ground state probability density at x."""
    return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-x**2 / (2 * sigma**2))

def qho_energy(sigma, m):
    """Compute QHO ground state energy E_0 = hbar^2 / (4 m sigma^2)."""
    return (HBAR**2) / (4 * m * sigma**2)

def validate_qho():
    """Validate Lambert W method for QHO with electron and trapped ion cases."""
    # Electron validation (Table S1, m = electron mass)
    print("QHO Validation (Electron):")
    print("Target σ (nm) | x_0 (nm) | y_0 (1/nm) | Calc. σ (nm) | E_0 (eV) | Rel. Error (%)")
    print("-" * 70)
    
    test_cases = [
        (0.4e-9, 0.1e-9),  # sigma = 0.4 nm, x_0 = 0.1 nm
        (0.8e-9, 0.5e-9),  # sigma = 0.8 nm, x_0 = 0.5 nm
        (1.6e-9, 1.0e-9),  # sigma = 1.6 nm, x_0 = 1.0 nm
    ]
    
    for sigma_true, x0 in test_cases:
        y0 = qho_pdf(x0, sigma_true)
        try:
            sigma_calc = compute_sigma(x0, y0, mu=0.0)
            e0_calc = qho_energy(sigma_calc, M_ELECTRON) / 1.60217662e-19  # J to eV
            rel_error = abs(sigma_calc - sigma_true) / sigma_true * 100
            print(f"{sigma_true*1e9:.1f} | {x0*1e9:.1f} | {y0*1e-9:.4f} | "
                  f"{sigma_calc*1e9:.4f} | {e0_calc:.5f} | {rel_error:.2f}")
        except ValueError as e:
            print(f"Error for sigma={sigma_true*1e9:.1f} nm: {e}")

    # Trapped ion (40Ca+) validation
    print("\nTrapped Ion (40Ca+) Validation:")
    m = M_CA_ION
    omega = 2 * np.pi * 1e6  # 1 MHz
    sigma_true = np.sqrt(HBAR / (2 * m * omega))  # ~10 nm
    x0 = 5e-9  # 5 nm
    y0 = qho_pdf(x0, sigma_true)  # ~1.59e7 1/m
    
    try:
        sigma_calc = compute_sigma(x0, y0, mu=0.0)
        e0_calc = qho_energy(sigma_calc, m) / 1.60217662e-19  # J to eV
        rel_error = abs(sigma_calc - sigma_true) / sigma_true * 100
        print(f"Target σ: {sigma_true*1e9:.1f} nm")
        print(f"x_0: {x0*1e9:.1f} nm, y_0: {y0:.2e} 1/m")
        print(f"Calculated σ: {sigma_calc*1e9:.1f} nm")
        print(f"E_0: {e0_calc:.2e} eV")
        print(f"Relative Error: {rel_error:.2f}%")
    except ValueError as e:
        print(f"Error for trapped ion: {e}")

if __name__ == "__main__":
    validate_qho()
```