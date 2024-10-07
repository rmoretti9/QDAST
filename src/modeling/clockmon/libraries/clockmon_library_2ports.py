import os
import pandas as pd
import re
import ast
from scipy import interpolate
from modeling.waveguides.libraries.waveguide_library import waveguide_library
import numpy as np
from qucat import Network,L,J,C,R
import numpy as np

def clockmon_library_2ports(deembed = 200, port_id = "0_2"):
    """
    Returns an interpolating function that maps coupler widths to capacitance matrices.
    Deembedding is performed using the waveguide library.
    
    Parameters
    ----------
    deembed : int, optional
        The deembedding distance in um. Default is 200 um.
    
    Returns
    -------
    library : scipy.interpolate.LinearNDInterpolator
        An interpolating function that maps coupler widths to capacitance matrices.
    """
    dir_name = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(dir_name, "clockmon_capacitance_library_" + port_id + "_sim_q3d_results.csv")
    df = pd.read_csv(file_path)
    cplr0_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[0])).values
    cplr2_widths = df['coupler_widths'].apply(lambda x: float(ast.literal_eval(x)[2])).values
    sweep_dim = len(cplr2_widths)
    CMatrix = np.zeros((sweep_dim, 4, 4))
    for i in range(sweep_dim):
        CMatrix[i] = df.iloc[i].values[2:].reshape(4,4)

    deembed = 200
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] = CMatrix[:, 0, 0] - wg_lib(deembed)
    CMatrix[:, 1, 1] = CMatrix[:, 1, 1] - wg_lib(deembed)

    library = interpolate.LinearNDInterpolator((cplr0_widths, cplr2_widths), CMatrix)
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
    network.append(J(3, 4, 'Lj'))
    cir = Network(network)
    Lj = 10e-9
    f, _, _, _ = cir.f_k_A_chi(Lj = Lj)
    c_sigma = 1/Lj/(2*np.pi*f[0])**2
    c_qr_1 = (CMatrix[0,2]*CMatrix[3,3] - CMatrix[0,3]*CMatrix[2,2]) / (CMatrix[2,2] + CMatrix[0,2] + CMatrix[3,3] + CMatrix[0,3])
    c_qr_2 = (CMatrix[1,2]*CMatrix[3,3] - CMatrix[1,3]*CMatrix[2,2]) / (CMatrix[2,2] + CMatrix[1,2] + CMatrix[3,3] + CMatrix[1,3])
    return c_sigma, c_qr_1, c_qr_2
