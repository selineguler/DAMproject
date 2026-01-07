from mace.calculators import mace_mp
from ase import build
from ase.units import GPa

def mace():
    calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')
    return calc
