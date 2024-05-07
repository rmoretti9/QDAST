import pandas as pd
import numpy as np
from scipy import interpolate
import re

def clockmon_library():
    df = pd.read_csv("clockmon_capacitance_library_sim_q3d_results.csv")
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
    library = interpolate.interp1d(coupler_widths, CMatrix, axis = 0)
    return library