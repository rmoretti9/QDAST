from math import pi
from qucat import Network,L,J,C,R
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re

from scipy.optimize import fsolve
from scqubits import Transmon

planck_h = 6.62607015e-34
e_charge = 1.60217663e-19
phi0 =  2.067833848e-15 # Wb

def get_Ic_from_Lj(Lj):
    return phi0/2/pi/Lj

def get_Lj_from_Ej(Ej):
    return (phi0/2/pi)**2 / Ej / planck_h

def get_Ej_from_Lj(Lj):
    return (phi0/2/pi)**2 / Lj / planck_h
