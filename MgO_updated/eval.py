import numpy as np
from eos import (
    murnaghan_energy,
    birch_murnaghan_energy_2nd,
    birch_murnaghan_energy_3rd,
)


def evaluate_eos(model, V, params):
    """
    Evaluate energy EOS at volumes V.
    params = dict containing required EOS params.
    """

    if model == "murnaghan":
        return murnaghan_energy(V, **params)

    if model == "bm2":
        return birch_murnaghan_energy_2nd(V, **params)

    if model == "bm3":
        return birch_murnaghan_energy_3rd(V, **params)

    raise ValueError(f"Unknown EOS {model}")

