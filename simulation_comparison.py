```python
"""
simulation_comparison.py

Compares the Lambert W method against least-squares fitting for estimating the
standard deviation (sigma) of a Gaussian PDF, as described in the PRL paper
"Single-Point Parameter Estimation for Gaussian Distributions Using the Lambert W
Function." Reproduces results from Supplemental Material (Section S-H) with 1%
Gaussian noise.

Usage:
    python simulation_comparison.py
"""

import numpy as np
from scipy.special import lambertw
from scipy.optimize import curve_fit
import time
from lambert_w_sigma import compute_sigma

def gaussian(x, sigma, A=1.0):
    """Gaussian PDF with mean mu=0."""
    return A / (sigma * np.sqrt(2 * np.pi)) * np.exp(-x**2 / (2 * sigma**2))

def run_lambert_w_simulation(x0, sigma_true, noise_level, n_trials=1000):
    """Run Lambert W method with noise, measure sigma and time."""
    y0_true = gaussian(x0, sigma_true)
    sigma_estimates = []
    start_time = time.time()
    
    for _ in range(n_trials):
        y0 = y0_true * (1 + np.random.normal(0, noise_level))
        try:
            sigma = compute_sigma(x0, y0, mu=0.0)
            sigma_estimates.append(sigma)
        except ValueError:
            continue
    
    end_time = time.time()
    mean_sigma = np.mean(sigma_estimates)
    std_sigma = np.std(sigma_estimates)
    comp_time = (end_time - start_time) / n_trials * 1e6  # us per trial
    
    return mean_sigma, std_sigma, comp_time

def run_least_squares_simulation(x_values, sigma_true, noise_level, n_trials=1000):
    """Run least-squares fitting with noise, measure sigma and time."""
    y_true = gaussian(x_values, sigma_true)
    sigma_estimates = []
    start_time = time.time()
    
    for _ in range(n_trials):
        y_noisy = y_true * (1 + np.random.normal(0, noise_level, size=len(x_values)))
        try:
            popt, _ = curve_fit(gaussian, x_values, y_noisy, p0=[sigma_true, 1.0])
            sigma_estimates.append(popt[0])
        except RuntimeError:
            continue
    
    end_time = time.time()
    mean_sigma = np.mean(sigma_estimates)
    std_sigma = np.std(sigma_estimates)
    comp_time = (end_time - start_time) / n_trials * 1e3  # ms per trial
    
    return mean_sigma, std_sigma, comp_time

def compare_methods():
    """Compare Lambert W and least-squares methods."""
    sigma_true = 1.0
    noise_level = 0.01  # 1% Gaussian noise
    n_trials = 1000
    
    # Lambert W: single-point
    x0 = 1.0  # Arbitrary point, away from mu=0
    lambert_mean, lambert_std, lambert_time = run_lambert_w_simulation(
        x0, sigma_true, noise_level, n_trials
    )
    
    # Least-squares: 10 points
    x_values = np.linspace(-3, 3, 10)
    ls_mean, ls_std, ls_time = run_least_squares_simulation(
        x_values, sigma_true, noise_level, n_trials
    )
    
    print("Simulation Comparison (sigma_true=1.0, 1% noise, 1000 trials):")
    print("-" * 60)
    print("Method         | Mean σ | Std. Dev. | Time per Trial")
    print("-" * 60)
    print(f"Lambert W      | {lambert_mean:.3f} | {lambert_std:.3f} | {lambert_time:.2f} us")
    print(f"Least-Squares  | {ls_mean:.3f} | {ls_std:.3f} | {ls_time:.2f} ms")

if __name__ == "__main__":
    compare_methods()
```