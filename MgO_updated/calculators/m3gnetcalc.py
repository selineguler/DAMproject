from ase import build
from ase.units import GPa
import matgl
from matgl.ext.ase import M3GNetCalculator

def m3gnet():
    matglpot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
    calc = M3GNetCalculator(potential=matglpot, stress_weight=GPa)
    return calc