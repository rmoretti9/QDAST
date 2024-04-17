from math import pi
import numpy as np
from scipy.special import ellipk

from math import sqrt, inf, tanh

epsilon_0 = 8.854e-12
mu_0 = 4 * np.pi * 1e-7
speed_of_light = 299792458  # m / s


def cpw_cl_ll(w, s, epsilon_r):

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
    c_l, l_l = C / l, L / l
    v_ph = 1 / np.sqrt(c_l * l_l)
    resonator_frequency = v_ph / 2 / l
    l_r = Z0 * pi / (4 * pi * resonator_frequency)
    c_r = 1 / ((2 * pi * resonator_frequency) ** 2 * l_r)
    return c_r, l_r


import numpy as np
from math import *
from scipy.special import ellipk

__all__ = ["kappa_in"]


def resonator_kappa(resonance_frequency, coupling_capacitance, z0, lambda_half):
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
    kappa = (1 / pi) * (omega_r**3.0) * (coupling_capacitance**2.0) * (z0**2.0)
    if not lambda_half:
        kappa = kappa * 2
    return kappa


def cpw_kinetic_inductance(
    trace_width: float, gap_width: float, film_thickness: float, pen_depth: float
):
    """Calculate kinetic inductance per unit length for a coplanar waveguide


    Based on formulas from:
    N. SIMONS  "Coplanar Waveguide Circuits, Components, and Systems" Section 12.4.3


    Args:
        trace_width: waveguide trace width (m)
        gap_width: gap width between trace and ground plane (m)
        film_thickness: thickness of the superconducting metal film (m)
        pen_depth: London penetration depth of the superconducting metal (m)


    Returns:
          Lk (float): kinetic inductance per unit length (H/m)
    """
    if pen_depth == 0.0:
        return 0.0
    S = trace_width
    W = gap_width
    k = S / (S + 2 * W)
    Kk = ellipk(k**2)
    t_pi = film_thickness / pi
    A = sqrt(4 * t_pi**2 + S**2) / 2 - t_pi
    B = S**2 / (4 * A)
    C = B - t_pi + sqrt(t_pi**2 + W**2)
    D = 2 * t_pi + C
    E = 1.7 / sinh(film_thickness / (2 * pen_depth))
    F = 0.4 / sqrt(((B / A) ** 2 - 1) * (1 - (B / D) ** 2))
    Lk = mu_0 * pen_depth * C * (E + F) / (4 * A * D * Kk)
    return Lk


def cpw_with_ground(
    eps_r: float,
    trace_width: float,
    gap_width: float,
    substrate_height: float = inf,
    chip_distance: float = inf,
    Lk: float = 0.0,
):
    """Calculate circuit parameters for a coplanar waveguide with bottom and top ground


    Based on analytic formula from:
    N. SIMONS  "Coplanar Waveguide Circuits, Components, and Systems" Section 3.2.1 (and 2.2.1)


    Args:
        eps_r: Relative dielectric constant of the substrate
        trace_width: waveguide trace width (m)
        gap_width: gap width between trace and ground plane (m)
        substrate_height: thickness of the substrate (m)
        chip_distance: distance between the bottom chip and top chip (m)
        Lk: kinetic inductance per unit length (H/m)


    Returns:
        (tuple): tuple containing:


        * c_eff (float): Effective speed of light in the waveguide (m/s)
        * eps_eff (float): Effective dielectric constant of the waveguide
        * Z0 (float): Characteristic impedance (Ohm)
        * Cs (float): Capacitance per unit length (F/m)
        * Ls (float): Inductance per unit length (H/m)
    """
    a = trace_width
    b = trace_width + 2 * gap_width

    if substrate_height == inf:
        k3 = a / b
    else:
        k3 = tanh(pi * a / (4.0 * substrate_height)) / tanh(
            pi * b / (4.0 * substrate_height)
        )

    if chip_distance == inf:
        k4 = a / b
    else:
        k4 = tanh(pi * a / (4.0 * chip_distance)) / tanh(pi * b / (4.0 * chip_distance))

    # Note: the definition of scipy.special.ellipk is different from the definition used in the reference
    # ellipk(k**2) here corresponds to K(k) in the reference.
    ratio_k3 = ellipk(k3**2) / ellipk(1 - k3**2)
    ratio_k4 = ellipk(k4**2) / ellipk(1 - k4**2)
    ratio_k1 = ratio_k3

    # Compute capacitance Cs and inductance Ls per unit length
    C_air = 2 * epsilon_0 * (ratio_k3 + ratio_k4)
    Cs = 2 * epsilon_0 * (eps_r - 1) * ratio_k1 + C_air
    Ls = 1.0 / (C_air * speed_of_light**2) + Lk  # total inductance (external + kinetic)

    # Compute characteristic impedance, effective speed of light, and effective permittivity using Cs and Ls.
    Z0 = sqrt(Ls / Cs)
    c_eff = 1.0 / sqrt(Ls * Cs)
    eps_eff = Ls * Cs * speed_of_light**2

    return c_eff, eps_eff, Z0, Cs, Ls


def cpw_kinetic_inductance(
    trace_width: float, gap_width: float, film_thickness: float, pen_depth: float
):
    """Calculate kinetic inductance per unit length for a coplanar waveguide


    Based on formulas from:
    N. SIMONS  "Coplanar Waveguide Circuits, Components, and Systems" Section 12.4.3


    Args:
        trace_width: waveguide trace width (m)
        gap_width: gap width between trace and ground plane (m)
        film_thickness: thickness of the superconducting metal film (m)
        pen_depth: London penetration depth of the superconducting metal (m)


    Returns:
          Lk (float): kinetic inductance per unit length (H/m)
    """
    if pen_depth == 0.0:
        return 0.0
    S = trace_width
    W = gap_width
    k = S / (S + 2 * W)
    Kk = ellipk(k**2)
    t_pi = film_thickness / pi
    A = sqrt(4 * t_pi**2 + S**2) / 2 - t_pi
    B = S**2 / (4 * A)
    C = B - t_pi + sqrt(t_pi**2 + W**2)
    D = 2 * t_pi + C
    E = 1.7 / sinh(film_thickness / (2 * pen_depth))
    F = 0.4 / sqrt(((B / A) ** 2 - 1) * (1 - (B / D) ** 2))
    Lk = mu_0 * pen_depth * C * (E + F) / (4 * A * D * Kk)
    return Lk


def penetration_depth(film_thickness: float, Tc: float):
    """Calculate London penetration depth from film thickness, sheet resistance, and critical temperature


    Based on BCS theory, Mattis-Bardeen formula, and Ginzburg-Landau penetration depth.
    For reference see https://arxiv.org/pdf/1408.4347.pdf or https://arxiv.org/pdf/1007.4187.pdf.


    Args:
        film_thickness: thickness of the superconducting metal film (m)
        Rsq: sheet resistance of the superconducting metal film (ohm / square)
        Tc: critical temperature of the superconducting metal (K)


    Returns:
          pen_depth (float): London penetration depth of the superconducting metal (m)
    """
    Rsq = 1.5e3
    return 1.05e-3 * sqrt(Rsq * film_thickness / Tc)
