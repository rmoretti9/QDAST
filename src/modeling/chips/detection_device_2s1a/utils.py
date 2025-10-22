"""Utility function module specific for chip 'detection_device_2s1a'."""

import pandas as pd
import skrf as rf
from modeling.utils.network import assemble_network
import numpy as np
from scipy.optimize import curve_fit

def get_transmission_line_tee_points(filename, tee_a, tee_b, qubit_name = None):

    if qubit_name == "s1":
        fl_traits = pd.read_csv(filename).values[0][1:]
    elif qubit_name == "s2":
        fl_traits = pd.read_csv(filename).values[1][1:]
    elif qubit_name == "a":
        fl_traits = pd.read_csv(filename).values[2][1:]
    else:
        ValueError("Specify a qubit name between 's1', 's2', 'a'.")

    transmission_line_tee_points = [
        float(fl_traits[0]), # lp to input capacitor
        fl_traits[1] + (tee_a + 2*tee_b)/2, # input capacitor to tee
        fl_traits[2] + fl_traits[3] + fl_traits[4] + (tee_a + 2*tee_b), # tee to split point
        fl_traits[5] + (tee_a + 2*tee_b)/2, # split point to lp
        30 + fl_traits[6] + fl_traits[7] # split point to ground

    ]
    return transmission_line_tee_points

class NetworkAnalysis:

    def __init__(self, options):
        self.options = options

    def s21_study(self):
        cnx = assemble_network(self.options)
        cir = rf.Circuit(cnx)
        ntw = cir.network
        return ntw
    
    def t1_purcell_readout(self, qb_frequency, c_sigma, with_purcell_filter = True):
        options = self.options
        if with_purcell_filter:
            options["type"] = "T1_Purcell_estimation"
        else:
            options["type"] = "T1_Purcell_estimation_simple"
        options["frequency"] = [qb_frequency - 0.01, qb_frequency + 0.01]
        options["n_points"] = 51
        cnx = assemble_network(self.options)
        cir = rf.Circuit(cnx)
        ntw = cir.network
        freq_span_s1 = np.linspace(options["frequency"][0], options["frequency"][1], options["n_points"])
        qb_idx = np.argmin(abs(freq_span_s1 - qb_frequency))
        Y_real = ntw.y[qb_idx, 0, 0].real
        t1_purcell_readout = c_sigma/Y_real
        return t1_purcell_readout

    def t1_purcell_coupler(self, qb_frequency, options_cplr):

        options_cplr_copy = options_cplr.copy()
        qubit_c = options_cplr_copy["qubit_c"][options_cplr_copy["qubit_to_estimate"]]

        cnx = assemble_network(options_cplr_copy)
        cir = rf.Circuit(cnx)

        ntw = cir.network
        freq_span_s1 = np.linspace(options_cplr_copy["frequency"][0], options_cplr_copy["frequency"][1], options_cplr_copy["n_points"])
        qb_idx = np.argmin(abs(freq_span_s1 - qb_frequency))
        Y_real = ntw.y[qb_idx, 0, 0].real
        t1_purcell_overall = qubit_c/Y_real
        
        options_cplr_copy["cc"] = [1e-22, 1e-22] # almost zero coupling to the bus coupler
        cnx = assemble_network(options_cplr_copy)
        cir = rf.Circuit(cnx)
        ntw = cir.network
        qb_idx = np.argmin(abs(freq_span_s1 - qb_frequency))
        Y_real = ntw.y[qb_idx, 0, 0].real
        t1_purcell_readout_only = qubit_c/abs(Y_real)

        t1_purcell_coupler_only = (1/t1_purcell_overall - 1/t1_purcell_readout_only)**(-1)
        return t1_purcell_coupler_only


def fit_resonance(freq, magnitude, Ql0, Qc0):
    def S21_model(f, f0, Ql, Qc, phi, a0, a1, a2, a3):
        """
        f    : frequency array
        f0   : resonator frequency
        Ql   : loaded Q
        Qc   : coupling Q (real, positive)
        phi  : rotation angle of complex S21
        a0-a3: cubic background coefficients (complex)
        """
        x = (f - f0)/f0
        # Lorentzian in complex form
        S21_res = 1 - (Ql/Qc) / (1 + 2j*Ql*x)
        # Apply rotation
        S21_res_rot = S21_res * np.exp(1j*phi)
        # Add cubic background
        bg = a0 + a1*x + a2*x**2 + a3*x**3
        return np.abs(S21_res_rot + bg)

    # Initial guess
    f0_guess = freq[np.argmin(magnitude)]
    phi_guess = 0.0
    a0_guess, a1_guess, a2_guess, a3_guess = 0,0,0,0
    p0 = [f0_guess, Ql0, Qc0, phi_guess, a0_guess, a1_guess, a2_guess, a3_guess]

    # Fit magnitude
    popt1, _ = curve_fit(S21_model, freq, magnitude, p0=p0)

    # Extract fitted parameters
    # f0_fit, Ql_fit, Qc_fit, phi_fit, a0_fit, a1_fit, a2_fit, a3_fit = popt1
    return popt1