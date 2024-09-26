import os
import pandas as pd
import re
import ast
from scipy import interpolate
from modeling.waveguides.libraries.waveguide_library import waveguide_library
import numpy as np

def clockmon_library_2ports(deembed = 200):
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
    file_path = os.path.join(dir_name, "clockmon_capacitance_library_0_2_sim_q3d_results.csv")
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
