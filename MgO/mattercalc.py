import torch
from loguru import logger
from ase.build import bulk
from ase.units import GPa
from mattersim.forcefield import MatterSimCalculator


def mattersim():
    calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth")
    return calc