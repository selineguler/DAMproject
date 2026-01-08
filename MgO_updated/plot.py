import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

from eos import (
    murnaghan_pressure,
    birch_murnaghan_pressure_2nd,
    birch_murnaghan_pressure_3rd,
    birch_murnaghan_pressure_4th,
    birch_murnaghan_pressure_5th,
    sjeos_pressure,
)


EOS_FAMILIES = {
    "murnaghan": murnaghan_pressure,
    "bm2": birch_murnaghan_pressure_2nd,
    "bm3": birch_murnaghan_pressure_3rd,
    "bm4": birch_murnaghan_pressure_4th,
    "bm5": birch_murnaghan_pressure_5th,
    "sjeos": sjeos_pressure,
}

BM_ONLY = ["bm2", "bm3", "bm4", "bm5"]

LINESTYLES = {
    "murnaghan": "-",
    "bm2": "--",
    "bm3": "-.",
    "bm4": ":",
    "bm5": (0, (3, 1, 1, 1)),
    "sjeos": (0, (5, 5)),
}

COLORS = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
]


############################################
# LOAD FILES
############################################

files = sys.argv[1:]
if len(files) == 0:
    raise ValueError("Provide CSV files as arguments.")


datasets = {}
for f in files:
    df = pd.read_csv(f)
    name = os.path.splitext(os.path.basename(f))[0]
    datasets[name] = df


############################################
# FIT EOS PARAMETERS PER MLIP
############################################

def get_eos_params(df, model_name):
    eos = df[(df["Type"] == "EOS") & (df["Model"] == model_name)]
    if len(eos) == 0:
        return None
    row = eos.iloc[0]
    return row


############################################
# ====== PLOT RAW ENERGY–VOLUME ===========
############################################

plt.figure(figsize=(6, 4))

for color, (name, df) in zip(COLORS, datasets.items()):
    raw = df[df["Type"] == "Raw"]
    plt.plot(raw["Volume"], raw["Energy"], "o-", label=name, color=color, markersize=3)

plt.xlabel("Volume (Å$^3$)")
plt.ylabel("Energy (eV)")
plt.title("MgO Energy–Volume Curve with ... Fit")
plt.legend()
plt.tight_layout()
plt.show()


############################################
# CREATE VOLUME GRID
############################################

all_vols = np.concatenate([
    df[df["Type"] == "Raw"]["Volume"].values
    for df in datasets.values()
])

vgrid = np.linspace(all_vols.min(), all_vols.max(), 400)


############################################
# ====== PLOT 1: MLIP vs EOS ===============
############################################

plt.figure(figsize=(7, 5))

for color, (name, df) in zip(COLORS, datasets.items()):

    # raw PV not plotted, only fits
    for eos_name, eos_func in EOS_FAMILIES.items():

        params = get_eos_params(df, eos_name if eos_name != "bm2" else "birchmurnaghan")

        if params is None:
            continue

        V0 = params["v0"]
        B0 = params["Bulk_Modulus_GPa"]
        B0p = 4.0  # if missing

        if eos_name == "murnaghan":
            P = eos_func(vgrid, V0, B0, 4.0)

        elif eos_name == "bm2":
            P = eos_func(vgrid, V0, B0)

        elif eos_name == "bm3":
            P = eos_func(vgrid, V0, B0, 4.0)

        elif eos_name == "bm4":
            P = eos_func(vgrid, V0, B0, 4.0, 0.0)

        elif eos_name == "bm5":
            P = eos_func(vgrid, V0, B0, 4.0, 0.0, 0.0)

        elif eos_name == "sjeos":
            # simple coefficients, placeholder
            P = eos_func(vgrid, 1, 1, 1, 0)

        else:
            continue

        plt.plot(
            vgrid,
            P,
            linestyle=LINESTYLES[eos_name],
            color=color,
            alpha=0.9,
            label=f"{name}-{eos_name}"
        )

plt.xlabel("Volume")
plt.ylabel("Pressure")
plt.title("Pressure–Volume: MLIPs + EOS Models")
plt.legend(fontsize=8, ncols=2)
plt.tight_layout()
plt.show()


############################################
# ====== PLOT 2: BM ORDER COMPARISON ======
############################################

plt.figure(figsize=(7, 5))

for color, (name, df) in zip(COLORS, datasets.items()):

    params = get_eos_params(df, "birchmurnaghan")
    if params is None:
        continue

    V0 = params["v0"]
    B0 = params["Bulk_Modulus_GPa"]

    for bm in BM_ONLY:
        func = EOS_FAMILIES[bm]

        if bm == "bm2":
            P = func(vgrid, V0, B0)

        elif bm == "bm3":
            P = func(vgrid, V0, B0, 4.0)

        elif bm == "bm4":
            P = func(vgrid, V0, B0, 4.0, 0.0)

        elif bm == "bm5":
            P = func(vgrid, V0, B0, 4.0, 0.0, 0.0)

        plt.plot(
            vgrid,
            P,
            linestyle=LINESTYLES[bm],
            color=color,
            label=f"{name}-{bm}"
        )

plt.xlabel("Volume")
plt.ylabel("Pressure")
plt.title("Pressure–Volume: Birch–Murnaghan Orders")
plt.legend(fontsize=8, ncols=2)
plt.tight_layout()
plt.show()

