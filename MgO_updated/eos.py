import numpy as np


# ===============================
# Helpers
# ===============================

def _eta(V, V0):
    """Finite strain variable."""
    return (V0 / V)**(2.0 / 3.0)


# ===============================
# Murnaghan EOS
# ===============================

def murnaghan_energy(V, V0, B0, B0p, E0):
    """
    Murnaghan EOS energy.
    """
    term = (B0 * V) / (B0p * (B0p - 1.0))
    return E0 + term * (((V0 / V)**B0p) * (B0p - 1.0) + 1.0) - (B0 * V0)/(B0p - 1.0)


def murnaghan_pressure(V, V0, B0, B0p):
    """
    Murnaghan EOS pressure.
    """
    return (B0 / B0p) * ((V0 / V)**B0p - 1.0)



# ===============================
# Birch–Murnaghan EOS (General)
# ===============================

def birch_murnaghan_energy_3rd(V, V0, B0, B0p, E0):
    """
    3rd-order Birch–Murnaghan EOS energy.
    """
    x = _eta(V, V0) - 1.0
    pref = (9.0 * V0 * B0) / 16.0
    return E0 + pref * (x**3 * B0p + x**2 * (6 - 4*_eta(V, V0)))


def birch_murnaghan_pressure_3rd(V, V0, B0, B0p):
    """
    3rd-order Birch–Murnaghan EOS pressure.
    """
    e = _eta(V, V0)
    return (3.0/2.0)*B0*(e**(7.0/2.0) - e**(5.0/2.0)) * \
           (1 + (3.0/4.0)*(B0p - 4.0)*(e - 1.0))



# ===============================
# 2nd-order Birch–Murnaghan
# ===============================

def birch_murnaghan_energy_2nd(V, V0, B0, E0):
    """
    2nd-order Birch–Murnaghan (B0' fixed = 4)
    """
    return birch_murnaghan_energy_3rd(V, V0, B0, B0p=4.0, E0=E0)


def birch_murnaghan_pressure_2nd(V, V0, B0):
    """
    2nd-order BM pressure.
    """
    return birch_murnaghan_pressure_3rd(V, V0, B0, B0p=4.0)



# ===============================
# 4th-order Birch–Murnaghan
# ===============================

def birch_murnaghan_energy_4th(V, V0, B0, B0p, B0pp, E0):
    """
    4th-order BM EOS energy.
    Based on Taylor expansion of finite strain.
    """
    f = 0.5 * (_eta(V, V0) - 1.0)

    return E0 + (9.0 * V0 * B0 / 2.0) * (
        f**2 +
        (B0p - 4.0)/3.0 * f**3 +
        ( (B0pp + B0p*(B0p-7.0) + 143.0/9.0) / 8.0 ) * f**4
    )


def birch_murnaghan_pressure_4th(V, V0, B0, B0p, B0pp):
    """
    4th-order BM pressure.
    """
    f = 0.5 * (_eta(V, V0) - 1.0)
    return 3*B0*f*(1+2*f)**(5.0/2.0) * (
        1 + (B0p-4)*f + (5/2)*(B0pp + B0p*(B0p-7) + 143/9)*f**2/3
    )



# ===============================
# 5th-order Birch–Murnaghan
# ===============================

def birch_murnaghan_energy_5th(V, V0, B0, B0p, B0pp, B0ppp, E0):
    """
    5th-order BM EOS energy.
    """
    f = 0.5 * (_eta(V, V0) - 1.0)

    return E0 + (9.0*V0*B0/2.0) * (
        f**2
        + (B0p-4)/3.0 * f**3
        + ( (B0pp + B0p*(B0p-7) + 143/9)/8.0 ) * f**4
        + (B0ppp/48.0)*f**5
    )


def birch_murnaghan_pressure_5th(V, V0, B0, B0p, B0pp, B0ppp):
    """
    5th-order BM pressure.
    """
    f = 0.5 * (_eta(V, V0) - 1.0)

    return 3*B0*f*(1+2*f)**(5.0/2.0) * (
        1
        + (B0p-4)*f
        + (5/2)*(B0pp + B0p*(B0p-7) + 143/9)/3.0 * f**2
        + (5/3)*B0ppp*f**3
    )



# ===============================
# SJEOS (Stabilized Jellium)
# ===============================

def sjeos_energy(V, a, b, c, d):
    """
    Stabilized Jellium EOS energy.
    Parameters are fit coefficients.
    """
    return (
        a/V**2 +
        b/V**(4/3) +
        c/V**(2/3) +
        d
    )


def sjeos_pressure(V, a, b, c, d):
    """
    SJEOS pressure = −dE/dV
    """
    return -(
        -2*a/V**3
        - (4/3)*b/V**(7/3)
        - (2/3)*c/V**(5/3)
    )



# ===============================
# Birch finite strain polynomial
# ===============================

def birch_energy(V, V0, a0, a1, a2, a3, E0):
    """
    Birch finite strain polynomial:
    E = sum a_n f^n
    """
    f = 0.5 * (_eta(V, V0) - 1.0)

    return E0 + (
        a0*f +
        a1*f**2 +
        a2*f**3 +
        a3*f**4
    )

