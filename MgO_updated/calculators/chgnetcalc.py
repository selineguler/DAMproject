from ase import build
from ase.units import GPa
import matgl
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import CHGNetCalculator

def chgnet():
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)
    return calc