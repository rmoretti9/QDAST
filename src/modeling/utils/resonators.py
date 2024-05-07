import skrf as rf
from math import pi
from skrf.media.cpw import CPW
import pandas as pd
import numpy as np
from scipy.special import ellipk

from math import sqrt, inf, tanh

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

def get_equivalent_lc(C, L, l, Z0):
    c_l, l_l = C/l, L/l
    v_ph = 1/np.sqrt(c_l*l_l)
    resonator_frequency = v_ph/2/l
    l_r = Z0*pi/(4*pi*resonator_frequency)
    c_r = 1/((2*pi*resonator_frequency)**2*l_r)
    return c_r, l_r

def loaded_quarter_wave_notch_resonator_kappa(frequency: float, Z0: float, ZL: float, C: float):
    """ Returns decay rate (kappa) of loaded quarter-wave resonator


    Derivation is based on approximating kappa = 2*frequency * |S12|^2 / |S11|^2 of following two-port system:
    ::


                          C
        port1----[Z0]----||----[ZL]----port2


    Args:
        frequency: Loaded resonator frequency (Hz)
        Z0: Characteristic impedance of the resonator line (Ohm)
        ZL: Characteristic impedance of the transmission line (Ohm).
        C: Coupling capacitance to transmission line (F)


    Returns:
        (float): External decay rate kappa (Hz), such that the FWHM line width in Hz is kappa/(2*pi)
    """
    return 32 * (C * pi) ** 2 * Z0 * ZL * frequency ** 3
