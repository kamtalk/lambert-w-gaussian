import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
import os

def error_amplification_factors(X):
    """Compute error amplification factors for given X = y0 * |x0 - mu| * sqrt(2 * pi * e)."""
    # Convert X to z = -2 * pi * y0^2 * (x0 - mu)^2
    z = -(X**2) / np.e  # Since X = y0 * |x0 - mu| * sqrt(2 * pi * e), z = -2 * pi * y0^2 * (x0 - mu)^2
    w0 = lambertw(z, k=0).real
    factor_y0 = np.abs(1 / (1 + w0))  # Amplification for delta y0 / y0
    factor_x0 = np.abs(w0 / (1 + w0))  # Amplification for delta x0 / (x0 - mu)
    return factor_y0, factor_x0

def plot_error_analysis():
    """Generate error amplification plot and save to output/error_analysis.pdf."""
    # Create output directory if it doesn't exist
    os.makedirs("output", exist_ok=True)
    
    # Generate X values (0 to just below 1, where z = -1/e)
    X = np.linspace(0.01, 0.999, 500)
    factor_y0, factor_x0 = error_amplification_factors(X)
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(X, factor_y0, 'b-', label=r'$|1/(1+W_0(z))|$')
    plt.plot(X, factor_x0, 'r-', label=r'$|W_0(z)/(1+W_0(z))|$')
    plt.xlabel(r'$X = y_0 |x_0 - \mu| \sqrt{2\pi e}$')
    plt.ylabel('Error Amplification Factor')
    plt.title('Error Amplification in Lambert $W$ Method')
    plt.grid(True)
    plt.legend()
    plt.yscale('log')  # Log scale to show divergence
    plt.tight_layout()
    
    # Save plot
    plt.savefig('output/error_analysis.pdf', format='pdf')
    plt.close()

if __name__ == "__main__":
    print("Generating error analysis plot...")
    plot_error_analysis()
    print("Saved to output/error_analysis.pdf")