import numpy as np


def log_likelihood_gaussian(
    V_ref,
    E_ref,
    V_model,
    E_model,
    sigma
):
    """
    Gaussian log-likelihood assuming

        E_model(V) = E_ref(V) + N(0, sigma^2)

    Inputs
    -------
    V_ref      : reference volume grid (1D array)
    E_ref      : reference energies evaluated on V_ref (1D array)
    V_model    : volume grid for model EOS (1D array)
    E_model    : model energies on V_model grid (1D array)
    sigma      : assumed Gaussian noise std (in eV)

    Returns
    -------
    logL       : float (Gaussian log-likelihood)
    """

    # Interpolate model energies onto reference grid
    spline = np.interp(V_ref, V_model, E_model)

    residuals = spline - E_ref

    chi2 = np.sum((residuals / sigma) ** 2)

    N = len(E_ref)

    logL = -0.5 * chi2 - N * np.log(sigma * np.sqrt(2.0 * np.pi))

    return logL

