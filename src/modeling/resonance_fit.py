from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

class ResFit():
    """Takes in input VNA scans containing a resonance in the form [freq, S21[dBm]] and returns the fit parameters."""
    def __init__(self, ntw, id_min, id_max, ideal = False, fit_type = "S21", poly_noise = True):
        self.xdata = ntw.f[id_min:id_max]
        if fit_type == "S21":
            self.ydata_dbm = ntw.s_db[id_min:id_max, 1, 0]
            self.ydata = abs(ntw.s[id_min:id_max, 1, 0])
        elif fit_type == "S11":
            self.ydata_dbm = ntw.s_db[id_min:id_max, 0, 0]
            self.ydata = abs(ntw.s[id_min:id_max, 0, 0])
        self.fmin= self.xdata[np.argmin(self.ydata)]
        self.xdata_centered = self.xdata - self.fmin
        self.ideal = ideal
        self.poly_noise = poly_noise
        self.fix_fmin()

    def fix_fmin(self):
        """fmin is a parameter that depends on the input data but must not vary during fit. fix_fmin() returns the fit function with fixed fmin."""
        if self.ideal:
            if self.poly_noise:
                # def fit_function(f, p0, p1, p3, p4, p5, p6, p7):
                #     fun = abs(p5*f + p6*f**2 + p7*f**3 +p0*(1 - np.exp(1j*p3)*p1*(p1**(-1)) / (1 + 2j * p1 * ( f - p4 ) / self.fmin )))
                #     return fun
                def fit_function(f, p0, p1, p3, p4, p5, p6, p7):
                    fun = abs(p5*f + p6*f**2 + p7*f**3 +p0*(1 - np.exp(1j*p3)/ (1 + 2j * p1 * ( f - p4 ) / self.fmin )))
                    return fun
            else:
                def fit_function(f, p0, p1, p3, p4, p5, p6, p7):
                    fun = abs(p0*(1 - np.exp(1j*p3)*p1*(p1**(-1)) / (1 + 2j * p1 * ( f - p4 ) / self.fmin )))
                    return fun
        else:
            if self.poly_noise:
                def fit_function(f, p0, p1, p2, p3, p4, p5, p6, p7):
                    fun = abs(p5*f + p6*f**2 + p7*f**3 +p0*(1 - np.exp(1j*p3)*p1*(p1**-1 - (p2)**-1) / (1 + 2j * p1 * ( f - p4 ) / self.fmin )))
                    # abs(p5*f + p6*f**2 + p7*f**3 +p0*(1 - np.exp(1j*p3)*p1*(p1**-1 - (p2)**-1) / (1 + 2j * p1 * ( f - p4 ) / self.fmin )))
                    return fun
            else:
                def fit_function(f, p0, p1, p2, p3, p4, p5, p6, p7):
                    fun = abs(p0*(1 - np.exp(1j*p3)*p1*(p1**-1 - (p2)**-1) / (1 + 2j * p1 * ( f - p4 ) / self.fmin )))
                    # abs(p5*f + p6*f**2 + p7*f**3 +p0*(1 - np.exp(1j*p3)*p1*(p1**-1 - (p2)**-1) / (1 + 2j * p1 * ( f - p4 ) / self.fmin )))
                    return fun
        self.fit_function = fit_function

    def fit_resonance(self, initial_params = None):
        self.popt, self.pcov = curve_fit(self.fit_function, self.xdata_centered, self.ydata, p0 = initial_params, maxfev = int(1e6))
        print(f"Resonant frequency: {self.popt[0] + self.fmin} Hz")
        print(f"Q total: {self.popt[1]}")
        if self.ideal:
            print(f"Q internal: set to inf")
            print(f"Complex phase: {self.popt[2]}")
        else:
            print(f"Q internal: {self.popt[2]}")
            print(f"Complex phase: {self.popt[3]}")
    
    def draw_resonance(self, initial_params):
        plt.plot(self.xdata * 1e-9, self.ydata_dbm)
        plt.plot(self.xdata * 1e-9, 20*np.log10(self.fit_function(self.xdata_centered, *initial_params)))
        plt.legend(["Data", "Fit"])
        plt.ylabel("S21 [dBm]")
        plt.xlabel("Frequency [GHz]")

    def plot_resonance(self, params = None):

        if params is None:
            params = self.popt
        plt.plot(self.xdata * 1e-9, self.ydata_dbm)
        plt.plot(self.xdata * 1e-9, 20*np.log10(self.fit_function(self.xdata_centered, *params)))
        plt.legend(["Data", "Fit"])
        plt.ylabel("S21 [dBm]")
        plt.xlabel("Frequency [GHz]")
