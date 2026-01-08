import numpy as np
import pandas as pd

volumes_m3gnet1 = np.arange(10.5, 16.0, 0.5)
volumes_m3gnet2 = np.arange(16.0, 21.7, 0.3)
volumes = np.concatenate((volumes_m3gnet1, volumes_m3gnet2))
volumes = np.append(volumes, 19.190)
volumes = np.sort(volumes)


# ---- reference EOS params (MgO-like) ----
V0  = 11.25     # Å^3 per formula unit
B0  = 160.0     # GPa  (human-readable)
B0p = 4.1
E0  = -320.00   # eV per formula unit


# ---- unit conversion ----
# 1 GPa = 6.241509e-3 eV / Å^3
GPA_TO_EV_A3 = 6.241509e-3
B0_eVA3 = B0 * GPA_TO_EV_A3


def bm3_energy(V, V0, B0_eVA3, B0p, E0):
    """
    Birch–Murnaghan 3rd-order EOS, returning energy in eV.
    V, V0 in Å^3 ; B0 in eV/Å^3
    """
    eta = (V0 / V) ** (2.0 / 3.0)
    x = eta - 1.0
    return (
        E0
        + (9.0 * V0 * B0_eVA3 / 16.0)
        * (x**3 * B0p + x**2 * (6.0 - 4.0 * eta))
    )


energies = bm3_energy(volumes, V0, B0_eVA3, B0p, E0)

df = pd.DataFrame({
    "Volume": volumes,
    "Energy": energies
})

df.to_csv("mgo.csv", index=False)

print("Wrote literature_MgO_eos.csv with", len(df), "rows")

