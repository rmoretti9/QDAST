import pandas as pd
import numpy as np
from scipy import interpolate
import re
import os

from modeling.waveguides.libraries.waveguide_library import waveguide_library

def clockmon_library(deembed = 200):
    """
    Returns an interpolating function that maps coupler widths to capacitance matrices.
    Deembedding is performed using the waveguide library.
    
    Parameters
    ----------
    deembed : int, optional
        The deembedding distance in um. Default is 200 um.
        
    Returns
    -------
    library : scipy.interpolate.interp1d
        An interpolating function that maps coupler widths to capacitance matrices.
    """
    dir_name = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(dir_name, "sweeps/clockmon_capacitance_library_sim_q3d_results.csv")
    df = pd.read_csv(file_path)
    coupler_widths = df["coupler_extent"].apply(lambda x: float(re.search(r'[-+]?\d*\.\d+|\d+', x).group()))
    CMatrix = df.iloc[:, 2:].values.reshape(-1, 3, 3)
    
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] -= wg_lib(deembed)
    library = interpolate.interp1d(coupler_widths, CMatrix, axis = 0)
    return library, coupler_widths, CMatrix


def clockmon_coupling_libraries(deembed = 200):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(dir_name, "sweeps/clockmon_capacitance_library_sim_q3d_results.csv")
    df = pd.read_csv(file_path)    
    coupler_widths_str = df["coupler_extent"].values

    coupler_widths = []
    pattern = r'[-+]?\d*\.\d+|\d+'

    for s in coupler_widths_str:
        match = re.search(pattern, s)
        if match:
            coupler_widths.append(float(match.group()))
        else:
            coupler_widths.append(None)
    sweep_dim = len(coupler_widths)
    CMatrix = np.zeros((sweep_dim, 3, 3))
    for i in range(sweep_dim):
        CMatrix[i] = df.iloc[i].values[2:].reshape(3,3)   
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] = CMatrix[:, 0, 0] - wg_lib(deembed)

    c_sigmas, c_qrs = get_csigma_cqr(CMatrix)
    
    coupler_width_given_c_qr = interpolate.interp1d(c_qrs, coupler_widths)
    c_sigma_given_coupler_width = interpolate.interp1d(coupler_widths, c_sigmas)
    return coupler_width_given_c_qr, c_sigma_given_coupler_width

def clockmon_cqr_to_ground(deembed = 200):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(dir_name, "sweeps/clockmon_capacitance_library_sim_q3d_results.csv")
    df = pd.read_csv(file_path)
    coupler_widths_str = df["coupler_extent"].values

    coupler_widths = []
    pattern = r'[-+]?\d*\.\d+|\d+'

    for s in coupler_widths_str:
        match = re.search(pattern, s)
        if match:
            coupler_widths.append(float(match.group()))
        else:
            coupler_widths.append(None)
    sweep_dim = len(coupler_widths)
    CMatrix = np.zeros((sweep_dim, 3, 3))
    for i in range(sweep_dim):
        CMatrix[i] = df.iloc[i].values[2:].reshape(3,3)
    
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] = CMatrix[:, 0, 0] - wg_lib(deembed)
    _, c_qr = get_csigma_cqr(CMatrix)
    
    library = interpolate.interp1d(c_qr, CMatrix[:, 0, 0], axis = 0)
    return library

def get_csigma_cqr(cmatrix):
    # C11 = cmatrix[0, 0]
    C12 = cmatrix[0, 1] if cmatrix.ndim == 2 else cmatrix[:, 0, 1]
    C13 = cmatrix[0, 2] if cmatrix.ndim == 2 else cmatrix[:, 0, 2]
    C22 = cmatrix[1, 1] if cmatrix.ndim == 2 else cmatrix[:, 1, 1]
    C23 = cmatrix[1, 2] if cmatrix.ndim == 2 else cmatrix[:, 1, 2]
    C33 = cmatrix[2, 2] if cmatrix.ndim == 2 else cmatrix[:, 2, 2]

    # Formulas adapted from https://qudev.phys.ethz.ch/static/content/science/Documents/semester/Burkhard_Simon_SemesterThesis_130211.pdf
    # Or similarly: https://www.nature.com/articles/s41534-020-0269-1.pdf
    c_sigma = ((C33 + C13)*(C22 + C12))/(C33 + C22 + C13 + C12) + C23
    beta = (C33*C12 - C22*C13)/((C33+C13)*(C22 + C12) + (C33 + C22 + C13 + C12)*C23)
    # Note that c_qr can also be written (perhaps more intuitively) as:
    # (C12*C33 - C13*C22) / (C22 + C12 + C33 + C13)
    c_qr = beta*c_sigma
    return c_sigma, c_qr