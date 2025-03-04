import os
import pandas as pd
import re
import ast
from scipy import interpolate
from modeling.waveguides.libraries.waveguide_library import waveguide_library
import numpy as np
from qucat import Network,L,J,C,R
import numpy as np

def clockmon_library_3ports(deembed = 200, port_id = "3_4_5"):
    """
    Returns an interpolating function that maps coupler widths to capacitance matrices.
    Deembedding is performed using the waveguide library.
    
    Parameters
    ----------
    deembed : int, optional
        The deembedding distance in um. Default is 200 um.
    port_id : str, optional
        The port id of the qubit couplers. Default is "3_4_5".

    Returns
    -------
    library : scipy.interpolate.LinearNDInterpolator
        An interpolating function that maps coupler widths to capacitance matrices.
    """
    dir_name = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(dir_name, "sweeps/clockmon_capacitance_library_" + port_id + "_sim_q3d_results.csv")
    df = pd.read_csv(file_path)
    cplr0_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[int(port_id.split("_")[0])])).values
    cplr1_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[int(port_id.split("_")[1])])).values
    cplr2_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[int(port_id.split("_")[2])])).values
    sweep_dim = len(cplr0_widths)
    CMatrix = np.zeros((sweep_dim, 5, 5))
    for i in range(sweep_dim):
        CMatrix[i] = df.iloc[i].values[2:].reshape(5,5)

    deembed = 200
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] = CMatrix[:, 0, 0] - wg_lib(deembed)
    CMatrix[:, 1, 1] = CMatrix[:, 1, 1] - wg_lib(deembed)
    CMatrix[:, 2, 2] = CMatrix[:, 2, 2] - wg_lib(deembed)

    library = interpolate.LinearNDInterpolator((cplr0_widths, cplr1_widths, cplr2_widths), CMatrix)
    return library

def clockmon_coupling_libraries(deembed = 200, port_id = "3_4_5"):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(dir_name, "sweeps/clockmon_capacitance_library_" + port_id + "_sim_q3d_results.csv")
    df = pd.read_csv(file_path)
    cplr0_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[int(port_id.split("_")[0])])).values
    cplr1_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[int(port_id.split("_")[1])])).values
    cplr2_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[int(port_id.split("_")[2])])).values
    sweep_dim = len(cplr0_widths)
    CMatrix = np.zeros((sweep_dim, 5, 5))
    for i in range(sweep_dim):
        CMatrix[i] = df.iloc[i].values[2:].reshape(5,5)
    c_qr_0, c_qr_1, c_qr_2 = get_cqr(CMatrix)

    points = np.column_stack((c_qr_0, c_qr_1, c_qr_2))
    values = np.column_stack((cplr0_widths, cplr1_widths, cplr2_widths))

    # Create the interpolator
    library = interpolate.LinearNDInterpolator(points, values)
    return library

    

def get_csigma_cqr(CMatrix):

    """
    Calculate C_sigma, C_qr_1 and C_qr_2 from the capacitance matrix of a 2-port clockmon qubit.
    
    Parameters
    ----------
    CMatrix : numpy array
        The capacitance matrix of the qubit.
    
    Returns
    -------
    c_sigma : float
        C_sigma in Farads.
    c_qr_1 : float
        C_qr_1 in Farads.
    c_qr_2 : float
        C_qr_2 in Farads.
    
    Notes
    -----
    C_sigma is the self-capacitance of the island and C_qr_1 and C_qr_2 are the effective capacitances between the coupler and the qubits. 
    The formula ignores the capacitance between couplers.
    """

    dim = CMatrix.shape[0]
    network = []
    for i in range(dim+ 1):
        for j in range(i+1, dim + 1):
            if i == 0:
                network.append(C(i, j, CMatrix[j-1, j-1]))
            else:
                network.append(C(i, j, CMatrix[i-1, j-1]))
    network.append(L(dim-1, dim, 'Lj'))
    cir = Network(network)
    Lj = 15e-9

    f = cir.eigenfrequencies(Lj = Lj)
    c_sigma = 1/Lj/(2*np.pi*max(f))**2
    c_qr_1, c_qr_2, c_qr_3 = get_cqr(CMatrix)
    return c_sigma, c_qr_1, c_qr_2, c_qr_3

def get_cqr(CMatrix):
    C14 = CMatrix[0, 3] if CMatrix.ndim == 2 else CMatrix[:, 0, 3]
    C15 = CMatrix[0, 4] if CMatrix.ndim == 2 else CMatrix[:, 0, 4]
    C24 = CMatrix[1, 3] if CMatrix.ndim == 2 else CMatrix[:, 1, 3]
    C25 = CMatrix[1, 4] if CMatrix.ndim == 2 else CMatrix[:, 1, 4]
    C34 = CMatrix[2, 3] if CMatrix.ndim == 2 else CMatrix[:, 2, 3]
    C35 = CMatrix[2, 4] if CMatrix.ndim == 2 else CMatrix[:, 2, 4]
    C44 = CMatrix[3, 3] if CMatrix.ndim == 2 else CMatrix[:, 3, 3]
    C55 = CMatrix[4, 4] if CMatrix.ndim == 2 else CMatrix[:, 4, 4]

    c_qr_1 = abs(C14*C55 - C15*C44) / (C44 + C14 + C55 + C15)
    c_qr_2 = abs(C24*C55 - C25*C44) / (C44 + C24 + C55 + C25)
    c_qr_3 = abs(C34*C55 - C35*C44) / (C44 + C34 + C55 + C35)
    return c_qr_1, c_qr_2, c_qr_3

def clockmon_cqr_to_ground(deembed = 200, port_id = "3_4_5", ground_id = 0):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(dir_name, "sweeps/clockmon_capacitance_library_" + port_id + "_sim_q3d_results.csv")
    df = pd.read_csv(file_path)
    cplr0_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[int(port_id.split("_")[0])])).values
    sweep_dim = len(cplr0_widths)
    CMatrix = np.zeros((sweep_dim, 5, 5))
    for i in range(sweep_dim):
        CMatrix[i] = df.iloc[i].values[2:].reshape(5,5)
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] = CMatrix[:, 0, 0] - wg_lib(deembed)
    CMatrix[:, 1, 1] = CMatrix[:, 1, 1] - wg_lib(deembed)
    CMatrix[:, 2, 2] = CMatrix[:, 2, 2] - wg_lib(deembed)

    c_qr_1, c_qr_2, c_qr_3 = get_cqr(CMatrix)
    library = interpolate.LinearNDInterpolator((c_qr_1, c_qr_2, c_qr_3), CMatrix[:, ground_id, ground_id])
    return library