from math import pi
from qucat import Network, L, J, C, R
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
from scqubits import Transmon
from modeling.utils.constants import planck_h, e_charge, phi0


def get_Ic_from_Lj(Lj):
    return phi0 / 2 / pi / Lj


def get_Lj_from_Ej(Ej):
    return (phi0 / 2 / pi) ** 2 / Ej / planck_h


def get_Ej_from_Lj(Lj):
    return (phi0 / 2 / pi) ** 2 / Lj / planck_h


def jaynes_cummings_g(EC, Ej, cqr, cr, rr_freq, lambda_half=False):
    """Calculate the coupling strength g between a transmon qubit and a resonator
    using the Jaynes-Cummings model.

    Parameters:
    EC (float): Charging energy of the transmon in Hz.
    Ej (float): Josephson energy of the transmon in Hz.
    cqr (float): Coupling capacitance between the transmon and the resonator in F.
    cr (float): Capacitance of the resonator in F.
    rr_freq (float): Resonator frequency in Hz.
    Returns:
    float: Coupling strength g in Hz.
    """
    g = (
        EC
        / e_charge
        * (Ej / EC / 2) ** (1 / 4)
        * cqr
        / cr
        * np.sqrt(2 * planck_h * rr_freq * cr)
    )
    if lambda_half:
        g = g / np.sqrt(2)
    return g
