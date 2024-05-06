import pandas as pd
import matplotlib.pyplot as plt
from qdast.modeling.utils.resonators import cpw_capacitance_inductance
from scipy import interpolate

def waveguide_library():
    df = pd.read_csv("waveguide_coplanar_library_sim_q3d_results.csv")
    C11 = df["C11"].values
    waveguide_length = df["wg_length"].values
    interp_function = interpolate.interp1d(waveguide_length, C11)
    return interp_function