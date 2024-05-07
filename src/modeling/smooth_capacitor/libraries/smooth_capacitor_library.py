import pandas as pd
import numpy as np
from scipy import interpolate

def make_library(filename):
    df = pd.read_csv(filename)
    finger_control = df["finger_control"].values
    sweep_dim = len(finger_control)
    CMatrix = np.zeros((sweep_dim, 2, 2))
    for i in range(sweep_dim):
        CMatrix[i] = df.iloc[i].values[2:].reshape(2,2)
    library = interpolate.interp1d(finger_control, CMatrix, axis = 0)
    return library

def digit_tee_library():
    filename = "digittee_capacitance_library_sim_q3d_results.csv"
    library = make_library(filename)
    return library

def smooth_capacitor_library():
    filename = "smooth_capacitor_capacitance_library_sim_output_results.csv"
    library = make_library(filename)
    return library