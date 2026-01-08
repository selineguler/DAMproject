import pandas as pd
import numpy as np

from likelihood import log_likelihood_gaussian
from eval import evaluate_eos


mlip_calculation = "out/chgnet.csv"
reference_file = "mgo.csv"

# NOTE: 0.02 eV (20 meV) is a reasonable starting value.
# You can tighten later if needed.
sigma = 2e-2   # eV

model = "bm3"  # here choose a model for EoS from evaluate_eos: currently bm2, bm3, murnaghan


df = pd.read_csv(mlip_calculation)

# get Birch–Murnaghan fit parameters
bm = df[(df["Type"] == "EOS") & (df["Model"] == "birchmurnaghan")].iloc[0]
GPA_TO_EV_A3 = 6.241509e-3
params = dict(
    V0=bm["v0"],
    B0=bm["Bulk_Modulus_GPa"] * GPA_TO_EV_A3, # this converts GPa to E/A3
    B0p=4.0,
    E0=bm["e0"],
)


# === define model volume grid (range around V0) ===
V_model = np.linspace(bm["v0"] * 0.9, bm["v0"] * 1.1, 400)

# evaluate model EOS energy
E_model = evaluate_eos(model, V_model, params)



# === load reference EOS data ===
df_ref = pd.read_csv(reference_file)

V_ref = df_ref["Volume"].values
E_ref = df_ref["Energy"].values


print("Before normalization:")
print("Model energy   min/max =", np.min(E_model), np.max(E_model))
print("Reference      min/max =", np.min(E_ref), np.max(E_ref))

# === NORMALIZE ENERGIES TO THEIR MINIMA ===
# This removes arbitrary energy offsets and ensures
# likelihood measures EOS shape differences only.
E_model = E_model - np.min(E_model)
E_ref   = E_ref   - np.min(E_ref)


print("After normalization:")
print("Model energy   min/max =", np.min(E_model), np.max(E_model))
print("Reference      min/max =", np.min(E_ref), np.max(E_ref))

# === compute log-likelihood ===
logL = log_likelihood_gaussian(
    V_ref,
    E_ref,
    V_model,
    E_model,
    sigma=sigma
)



print(f"MLIP file: {mlip_calculation}")
print(f"Reference: {reference_file}")
print(f"EOS model: {model}")
print(f"Assumed sigma: {sigma} eV")
print(f"Log-likelihood: {logL:.6f}")

