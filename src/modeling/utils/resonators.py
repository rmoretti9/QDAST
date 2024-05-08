import skrf as rf
from math import pi
from skrf.media.cpw import CPW
import pandas as pd
import numpy as np
from scipy.special import ellipk

from math import sqrt, inf, tanh

def cpw_cl_ll(w, s, epsilon_r):
    epsilon_0 = 8.854e-12
    
    # Permeability of free space
    mu_0 = 4 * np.pi * 1e-7
    
    # Effective dielectric constant
    epsilon_eff = (epsilon_r + 1) / 2
    
    # Calculation of the elliptic integrals
    k = w / (w + 2 * s)
    k_prime = np.sqrt(1 - k**2)
    
    K_k = ellipk(k)
    K_k_prime = ellipk(k_prime)
    
    # Capacitance per unit length
    C_l = 4 * epsilon_0 * epsilon_eff * K_k / K_k_prime
    
    # Inductance per unit length
    L_l = mu_0 * K_k_prime / (4 * K_k)
    return C_l, L_l

def cpw_capacitance_inductance(w, s, l, epsilon_r):
    """
    Calculates the capacitance and inductance of a CPW resonator.
    
    Parameters:
    w (float): Width of the center conductor (m)
    s (float): Gap between the center conductor and ground planes (m)
    l (float): Length of the CPW resonator (m)
    epsilon_r (float): Relative permittivity of the substrate
    
    Returns:
    C (float): Capacitance of the CPW resonator (F)
    L (float): Inductance of the CPW resonator (H)
    """
    
    C_l, L_l = cpw_cl_ll(w, s, epsilon_r)
    C = C_l * l
    L = L_l * l
    
    return C, L

def get_equivalent_lc(C, L, l, Z0):
    c_l, l_l = C/l, L/l
    v_ph = 1/np.sqrt(c_l*l_l)
    resonator_frequency = v_ph/2/l
    l_r = Z0*pi/(4*pi*resonator_frequency)
    c_r = 1/((2*pi*resonator_frequency)**2*l_r)
    return c_r, l_r

import numpy as np
from math import *
from scipy.special import ellipk

__all__ = ['kappa_in']


def resonator_kappa(resonance_frequency, coupling_capacitance, z0):
    """Key References:

    D. Schuster, Ph.D. Thesis, Yale University (2007)
    https://rsl.yale.edu/sites/default/files/files/RSL_Theses/SchusterThesis.pdf

    T. McConkey, Ph.D. Thesis, University of Waterloo (2018)
    https://uwspace.uwaterloo.ca/bitstream/handle/10012/13464/McConkey_Thomas.pdf?sequence=3&isAllowed=y

    Mohebbi and Majedi, Superconducting Science and Technology 22, 125028 (2009)
    https://iopscience.iop.org/article/10.1088/0953-2048/22/12/125028/meta

    P. Krantz, et al. Physical Review Applied 6, 021318 (2019)
    https://aip.scitation.org/doi/10.1063/1.5089550
    """

    # Calculation of kappa
    omega_r = resonance_frequency * 2 * pi
    kappa = (2 / pi) * (omega_r**3.0) * (coupling_capacitance**2.0) * (z0**2.0)
    return kappa