import ase
from ase.build import bulk

from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
def orb():
    device="cpu" 
    orbff = pretrained.orb_v3_conservative_inf_omat(
    device=device,
    precision="float32-high",  
    )
    calc = ORBCalculator(orbff, device=device)
    return calc