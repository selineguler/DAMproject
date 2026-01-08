import os
import sys
import numpy as np
import pandas as pd
from ase.build import bulk
from ase.eos import EquationOfState
from ase.units import kJ
from volumes import volumes_m3gnet
from scipy.interpolate import InterpolatedUnivariateSpline

from eos import (
    murnaghan_pressure,
    birch_murnaghan_pressure_2nd,
    birch_murnaghan_pressure_3rd,
    birch_murnaghan_pressure_4th,
    birch_murnaghan_pressure_5th,
    sjeos_pressure,
)

def get_calculator(conda):
    if conda == "m3gnet":
        from calculators.m3gnetcalc import m3gnet
        return m3gnet()

    if conda == "mace":
        from calculators.macecalc import mace
        return mace()

    if conda == "sevenn":
        from calculators.sevenncalc import sevenn
        return sevenn()

    if conda == "mattersim":
        from calculators.mattercalc import mattersim
        return mattersim()

    if conda == "orb":
        from calculators.orbcalc import orb
        return orb()

    if conda == "chgnet":
        from calculators.chgnetcalc import chgnet
        return chgnet()


def run_eos(fout, conda):

    calculator = get_calculator(conda)
    volumes = volumes_m3gnet
    vfine = np.linspace(volumes.min(), volumes.max(), 500) # why we define Vfine?

    rows = []
    energies = []

    for vol in volumes:
        stru = bulk('MgO', crystalstructure='rocksalt', a=(4*vol)**(1/3))
        stru.calc = calculator

        E = stru.get_potential_energy()
        energies.append(E)

        rows.append({
            "Type": "Raw",
            "Volume": vol,
            "Energy": E,
        })

    df = pd.DataFrame(rows)

    spline = InterpolatedUnivariateSpline(volumes, energies, k=3)
    Einterp = spline(vfine) # why we define interpolated?

    interp_df = pd.DataFrame({
        "Type": "Interpolated",
        "Volume": vfine,
        "Energy": Einterp
    })


    eos_models = [
        ("birchmurnaghan",),
        ("murnaghan",),
        ("sjeos",),
        ("taylor",),
        ("birch",),
    ]

    eos_rows = []

    for (name,) in eos_models:
        eos = EquationOfState(volumes, energies, eos=name)
        v0, e0, B = eos.fit()

        eos_rows.append({
            "Type": "EOS",
            "Model": name,
            "v0": v0,
            "e0": e0,
            "Bulk_Modulus_GPa": B / kJ * 1e24
        })

    eos_df = pd.DataFrame(eos_rows)

    pv_rows = []

    # helper: get v0 and B0 from BM fit
    bm_row = eos_df[eos_df["Model"] == "birchmurnaghan"].iloc[0]
    v0 = bm_row["v0"]
    B0 = bm_row["Bulk_Modulus_GPa"]

    # Assume defaults for higher derivatives
    B0p = 4.0
    B0pp = 0.0
    B0ppp = 0.0

    for V in vfine:

        # --- Murnaghan ---
        P = murnaghan_pressure(V, v0, B0, B0p)
        pv_rows.append({
            "Type": "PV",
            "Model": "murnaghan",
            "Volume": V,
            "Pressure": P,
        })

        # --- BM2 ---
        P = birch_murnaghan_pressure_2nd(V, v0, B0)
        pv_rows.append({
            "Type": "PV",
            "Model": "bm2",
            "Volume": V,
            "Pressure": P,
        })

        # --- BM3 ---
        P = birch_murnaghan_pressure_3rd(V, v0, B0, B0p)
        pv_rows.append({
            "Type": "PV",
            "Model": "bm3",
            "Volume": V,
            "Pressure": P,
        })

        # --- BM4 ---
        P = birch_murnaghan_pressure_4th(V, v0, B0, B0p, B0pp)
        pv_rows.append({
            "Type": "PV",
            "Model": "bm4",
            "Volume": V,
            "Pressure": P,
        })

        # --- BM5 ---
        P = birch_murnaghan_pressure_5th(V, v0, B0, B0p, B0pp, B0ppp)
        pv_rows.append({
            "Type": "PV",
            "Model": "bm5",
            "Volume": V,
            "Pressure": P,
        })

        # --- SJEOS ---
        # NOTE: just placeholder coeffs unless fitted separately
        P = sjeos_pressure(V, 1.0, 1.0, 1.0, 0.0)
        pv_rows.append({
            "Type": "PV",
            "Model": "sjeos",
            "Volume": V,
            "Pressure": P,
        })

    pv_df = pd.DataFrame(pv_rows)

    # ===========================
    # Combine all
    # ===========================
    out = pd.concat(
        [df, eos_df, interp_df, pv_df],
        ignore_index=True
    )

    out.to_csv(os.path.join(fout, f"{conda}.csv"), index=False)


if __name__ == "__main__":
    fout = sys.argv[1]
    conda = sys.argv[2]
    run_eos(fout, conda)

