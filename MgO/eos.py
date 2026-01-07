
import os
import sys
from scipy.constants import angstrom            ## one Angstrom in meters
from scipy.constants import atmosphere          ## standard atmosphere in pascals
from scipy.constants import Boltzmann           ## Boltzmann constant
from scipy.constants import electron_mass       ## electron mass
from scipy.constants import electron_volt       ## one electron volt in Joules
from scipy.constants import giga
from scipy.constants import hbar                ## h/2pi
from scipy.constants import nano
from scipy.constants import milli
from scipy.constants import physical_constants
from scipy.constants import Planck              ## the Planck constant h
from scipy.optimize import curve_fit
import sympy as sp
from ase import Atom, Atoms, build
from ase.build import stack, surface
from ase.cell import Cell
from ase.constraints import FixAtoms, FixSymmetry
from ase.dft.kpoints import BandPath
from ase.filters import FrechetCellFilter
from ase.formula import Formula
from ase.geometry import get_distances
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.lattice import HEX
from ase.lattice.cubic import FaceCenteredCubic
from ase.optimize import FIRE
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.spacegroup.symmetrize import check_symmetry, refine_symmetry
from ase.units import GPa
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.eos import calculate_eos
from ase.units import kJ
from ase.visualize import view
import datetime
from importlib.metadata import version, PackageNotFoundError
import importlib.util
from io import StringIO
from ase.eos import EquationOfState
import numpy as np
from scipy.interpolate import bisplev, InterpolatedUnivariateSpline, RectBivariateSpline
from scipy.optimize import minimize, minimize_scalar
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
import pandas as pd
from volumes import volumes_m3gnet, volumes_grace

fout = sys.argv[1]
conda = sys.argv[2]


if conda == "m3gnet":
    from calculators import m3gnet
    calculator = m3gnet()

if conda == "chgnet":
    from calculators import chgnet
    calculator = chgnet()

if conda == "grace":
    from calculators import grace
    calculator = grace()

if conda == "mace":
    from calculators import mace
    calculator = mace()

if conda == "sevenn":
    from sevenncalculator import sevenn
    calculator = sevenn()

if conda == "mattersim":
    from mattercalc import mattersim
    calculator = mattersim()

volumes = volumes_m3gnet
vfine = np.linspace(volumes.min(), volumes.max(), 500)

rows = []
energies = []

for vol in volumes:
    stru = bulk('MgO', crystalstructure='rocksalt', a=np.power(4*vol, 1/3.0))
    stru.calc = calculator

    E = stru.get_potential_energy()
    energies.append(E)

    cell = np.array(stru.cell).reshape(-1)

    rows.append({
        "Type": "Raw",
        "Volume": vol,
        "Energy": E,
        "Periodic Boundary Conditions": stru.pbc,
        "a_x": cell[0], "a_y": cell[1], "a_z": cell[2],
        "b_x": cell[3], "b_y": cell[4], "b_z": cell[5],
        "c_x": cell[6], "c_y": cell[7], "c_z": cell[8],
    })

df = pd.DataFrame(rows)

spline = InterpolatedUnivariateSpline(volumes, energies, k=3)
Einterpolated = spline(vfine) 

eos_1     = EquationOfState(volumes, energies, eos='birchmurnaghan')
eos_2    = EquationOfState(volumes, energies, eos='murnaghan')
eos_3    = EquationOfState(volumes, energies, eos='sjeos')
eos_4    = EquationOfState(volumes, energies, eos='taylor')
eos_5    = EquationOfState(volumes, energies, eos='birch')

v0_1, e0_1, B_1 = eos_1.fit()
v0_2, e0_2, B_2 = eos_2.fit()   
v0_3, e0_3, B_3 = eos_3.fit()
v0_4, e0_4, B_4 = eos_4.fit()
v0_5, e0_5, B_5 = eos_5.fit()

eos_summary = pd.DataFrame([
    {"Type": "EOS", "Model": "birchmurnaghan", "v0": v0_1, "e0": e0_1, "Bulk_Modulus_GPa": B_1 / kJ * 1e24},
    {"Type": "EOS", "Model": "murnaghan",      "v0": v0_2, "e0": e0_2, "Bulk_Modulus_GPa": B_2 / kJ * 1e24},
    {"Type": "EOS", "Model": "sjeos",          "v0": v0_3, "e0": e0_3, "Bulk_Modulus_GPa": B_3 / kJ * 1e24},
    {"Type": "EOS", "Model": "taylor",         "v0": v0_4, "e0": e0_4, "Bulk_Modulus_GPa": B_4 / kJ * 1e24},
    {"Type": "EOS", "Model": "birch",          "v0": v0_5, "e0": e0_5, "Bulk_Modulus_GPa": B_5 / kJ * 1e24},
])

#eos_summary.loc[eos_summary["Model"] == "birchmurnaghan", "Bulk_Prime"] = eos_1.eos_parameters[2]
#eos_summary.loc[eos_summary["Model"] == "murnaghan",      "Bulk_Prime"] = eos_2.eos_parameters[2]
#eos_summary.loc[eos_summary["Model"] == "sjeos",          "Bulk_Prime"] = eos_3.eos_parameters[2]
#eos_summary.loc[eos_summary["Model"] == "taylor",         "Bulk_Prime"] = eos_4.eos_parameters[2]
#eos_summary.loc[eos_summary["Model"] == "birch",          "Bulk_Prime"] = eos_5.eos_parameters[2]

interpolated_df = pd.DataFrame({
    "Type": "Interpolated",
    "Volume": vfine,
    "Energy": Einterpolated,
})

out_df = pd.concat([df, eos_summary, interpolated_df], ignore_index=True)

out_df.to_csv(os.path.join(fout, f"{conda}.csv"), index=False)

