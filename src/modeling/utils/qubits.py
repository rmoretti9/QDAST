from math import pi
from modeling.clockmon.libraries.clockmon_library import clockmon_library
from qucat import Network,L,J,C,R
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skrf.media.cpw import CPW
import re

from scipy.optimize import fsolve
from scqubits import Transmon

planck_h = 6.62607015e-34
e_charge = 1.60217663e-19
phi0 =  2.067833848e-15 # Wb

def get_csigma_cqr(cmatrix):
    # C11 = cmatrix[0, 0]
    C12 = cmatrix[0, 1] if cmatrix.ndim == 2 else cmatrix[:, 0, 1]
    C13 = cmatrix[0, 2] if cmatrix.ndim == 2 else cmatrix[:, 0, 2]
    C22 = cmatrix[1, 1] if cmatrix.ndim == 2 else cmatrix[:, 1, 1]
    C23 = cmatrix[1, 2] if cmatrix.ndim == 2 else cmatrix[:, 1, 2]
    C33 = cmatrix[2, 2] if cmatrix.ndim == 2 else cmatrix[:, 2, 2]

    # Formulas adapted from https://qudev.phys.ethz.ch/static/content/science/Documents/semester/Burkhard_Simon_SemesterThesis_130211.pdf
    c_sigma = ((C33 + C13)*(C22 + C12))/(C33 + C22 + C13 + C12) + C23
    beta = (C33*C12 - C22*C13)/((C33+C13)*(C22 + C12) + (C33 + C22 + C13 + C12)*C23)
    # Note that c_qr can also be written (perhaps more intuitively) as:
    # (C12*C33 - C13*C22) / (C22 + C12 + C33 + C13)
    c_qr = beta*c_sigma
    return c_sigma, c_qr


def get_Ic_from_Lj(Lj):
    return phi0/2/pi/Lj

def get_Lj_from_Ej(Ej):
    return (phi0/2/pi)**2 / Ej / planck_h

def get_Ej_from_Lj(Lj):
    return (phi0/2/pi)**2 / Lj / planck_h