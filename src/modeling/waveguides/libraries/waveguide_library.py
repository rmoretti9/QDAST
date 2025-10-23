import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import os


def waveguide_library():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(dir_name, "waveguide_coplanar_library_sim_q3d_results.csv")
    df = pd.read_csv(file_path)
    C11 = df["C11"].values
    waveguide_length = df["wg_length"].values
    interp_function = interpolate.interp1d(waveguide_length, C11)
    return interp_function
