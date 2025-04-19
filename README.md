Lambert W Gaussian Parameter Estimation
This private repository contains Python code to reproduce results for the Physical Review Letters paper, "Single-Point Parameter Estimation for Gaussian Distributions Using the Lambert ( W ) Function." The method derives the standard deviation (\sigma) of a normalized Gaussian probability density function (GPDF) from a single measurement ((x_0, y_0)), given the mean (\mu), using the Lambert ( W ) function (Eq. 1 in the paper). It enables microsecond-scale parameter extraction for real-time and data-scarce experiments in quantum physics, gravitational-wave astronomy, and more.
Project Structure

lambert_w_sigma.py: Computes (\sigma) using the Lambert ( W ) function (Eq. 1).
qho_validation.py: Validates the method for quantum harmonic oscillator (QHO) ground state, recovering (\sigma) and energy (E_0) (Table I).
maxwell_boltzmann.py: Validates the method for 1D Maxwell-Boltzmann velocity distribution, estimating temperature (T) (Table I).
simulation_comparison.py: Compares Lambert ( W ) method performance against least-squares fitting, highlighting speed and accuracy.
error_analysis.py: Computes error propagation (Eq. 3) and generates output/error_analysis.pdf to visualize noise sensitivity.
requirements.txt: Python dependencies for running the scripts.
output/: Directory containing error_analysis.pdf.

Setup

Clone the repository (private access required):
git clone <private-repo-url>
cd lambert-physics


Install dependencies:
pip install -r requirements.txt



Usage

Run QHO validation:
python qho_validation.py


Run Maxwell-Boltzmann validation:
python maxwell_boltzmann.py


Compare methods:
python simulation_comparison.py


Generate error analysis:
python error_analysis.py

Output: output/error_analysis.pdf


Experimental Data
Results are validated using experimental data:

Trapped ion ((^{40}\text{Ca}^+)): QHO ground state, recovering (\sigma \approx 10 , \text{nm}).
Ethanol 1H-NMR (SDBS No. 1300): CH₃ peak for Gaussian fitting, processed via lambert_w_sigma.py.

Notes

Requires Python 3.8+.
Contact a.researcher@institution.edu for access or inquiries.
Refer to the PRL paper for theoretical details and derivations.

Acknowledgments
Supported by [Funding Source, optional]. Computations use SciPy [Vir20].
