import skrf as rf
from math import pi
from skrf.media.cpw import CPW
import pandas as pd
import numpy as np
from scipy.special import ellipk

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
    
    # Permittivity of free space
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
    
    # Total capacitance and inductance
    C = C_l * l
    L = L_l * l
    
    return C, L