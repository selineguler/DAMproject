from ase.units import GPa
from fairchem.core.calculators.ocp import OCPCalculator



def grace():
    calc = OCPCalculator(checkpoint_path="eqV2_31M_omat_mp_salex.pt")
    return calc