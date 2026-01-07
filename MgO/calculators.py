from mace.calculators import mace_mp
from ase import build
from ase.units import GPa
import matgl
from matgl.ext.ase import M3GNetCalculator
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import CHGNetCalculator
from sevenn.calculator import SevenNetCalculator

def sevenn():    
    calc = SevenNetCalculator(
        model="/path/to/7net-omni",
        modal='mpa',
        enable_cueq=False,
        enable_flash=False
    )
    return calc

def mace():
    calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')
    return calc

def m3gnet():
    matglpot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
    calc = M3GNetCalculator(potential=matglpot, stress_weight=GPa)
    return calc

def chgnet():
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)
    return calc