import pandas as pd
import numpy as np
from scipy import interpolate
import os
from modeling.waveguides.libraries.waveguide_library import waveguide_library

dir_name = os.path.dirname(os.path.abspath(__file__))

def make_library(filename, deembed = [300, 200]):
    df = pd.read_csv(filename)
    finger_control = df["finger_control"].values
    CMatrix = get_cmatrices(filename, deembed)
    library = interpolate.interp1d(finger_control, CMatrix, axis = 0)
    return library

def digit_tee_library(deembed = [300, 200]):
    filename = os.path.join(dir_name, "digittee_capacitance_library_sim_q3d_results.csv")
    library = make_library(filename, deembed)
    return library

def digit_tee_ck_to_ground(deembed = [300, 200]):
    filename = os.path.join(dir_name, "digittee_capacitance_library_sim_q3d_results.csv")
    CMatrix = get_cmatrices(filename, deembed)
    library = interpolate.interp1d(CMatrix[:, 0, 1], CMatrix[:, 1, 1], axis = 0)
    return library

def get_ck(deembed = [300, 200], type = "digittee"):
    if type == "digittee":
        filename = os.path.join(dir_name, "digittee_capacitance_library_sim_q3d_results.csv")
    elif type == "smooth_capacitor":
        filename = os.path.join(dir_name, "smooth_capacitor_capacitance_library_sim_output_results.csv")
    df = pd.read_csv(filename)
    finger_control = df["finger_control"].values
    CMatrix = get_cmatrices(filename, deembed)
    library = interpolate.interp1d(CMatrix[:, 0, 1], finger_control, axis = 0)
    return library    

def get_cmatrices(filename, deembed = [300, 200]):
    df = pd.read_csv(filename)
    finger_control = df["finger_control"].values
    sweep_dim = len(finger_control)
    CMatrix = np.zeros((sweep_dim, 2, 2))
    for i in range(sweep_dim):
        CMatrix[i] = df.iloc[i].values[2:].reshape(2,2)

    wg_lib = waveguide_library()
    if isinstance(deembed, list):
        CMatrix[:, 0, 0] = CMatrix[:, 0, 0] - wg_lib(deembed[0])
        CMatrix[:, 1, 1] = CMatrix[:, 1, 1] - wg_lib(deembed[1])
    else:
        CMatrix[:, 0, 0] = CMatrix[:, 0, 0] - wg_lib(deembed)
        CMatrix[:, 1, 1] = CMatrix[:, 1, 1] - wg_lib(deembed)
    return CMatrix

def smooth_capacitor_library(deembed = 200):
    filename = os.path.join(dir_name, "smooth_capacitor_capacitance_library_sim_output_results.csv")
    library = make_library(filename, deembed)
    return library