"""Various calculations to estimate cross-resonance and spurious terms, 
assuming two transmon qubits coupled through a cavity."""
import numpy as np
import scqubits as scq
from scqubits import HilbertSpace
from scipy.optimize import fsolve

planck_h = 6.62607015e-34
e_charge = 1.60217663e-19
phi0 =  2.067833848e-15

class CRModel():
    """Cross-resonance model."""
    def __init__(self, w1, w2, wc, a1, a2, g1, g2, C1, C2, use_dressed_freq, mode_a):
        self.w1 = w1
        self.w2 = w2
        self.wc = wc
        self.a1_forced = a1
        self.a2_forced = a2
        self.g1 = g1
        self.g2 = g2
        self.C1 = C1
        self.C2 = C2
        self.diagonalize()
        self.J = self.compute_J(dressed = use_dressed_freq)

    def compute_J(self, dressed = False):
        if dressed:
            w1 = self.w1_dressed
            w2 = self.w2_dressed
        else:
            w1 = self.w1
            w2 = self.w2

        numerator = self.g2 * self.g1 * (w1 + w2 - 2 * self.wc)
        denominator = 2 * (w1 - self.wc) * (w2 - self.wc)
        J = numerator / denominator
        self.J = J
    
    def zz_strength(self, drive_strength, mode_a = "numeric"):
        D = self.w1 - self.w2
        if mode_a == "numeric":
            a1, a2 = self.a1, self.a2
        elif mode_a == "forced":
            a1, a2 = self.a1_forced, self.a2_forced
        elif mode_a == "dressed":
            a1, a2 = self.a1_dressed, self.a2_dressed

        self.zz_static = self.J**2 * (a1 + a2)/ (a1 + D) / (D - a2)
        self.zz_driven = self.J**2 / 2 / (a1 + D)**2 * drive_strength**2* (
            (a1**3 - 2*a1*D**2 - 2*D**3) / (a1*D**2 * (a2 -D)) + 
            0.5* ( 
                4*(3*a1 + D)*(a1**2 + a1*D + D**2) / (D**2 * (2*a1 + D)**2) -
                16*D / (3*a1**2 + 8*a1*D + 4*D**2)
                ) + 
            2*a1 / (D*a2) - 2*(a1 + D) / (a1 + D - a2)**2 - 2*(a1 + D) / a1 / (a1 + D - a2)
            )
        self.zz = self.zz_static + self.zz_driven
    
    def zx_strength(self, drive_strength):
        D = self.w1 - self.w2
        a1, a2 = self.a1, self.a2
        self.zx = -self.J*drive_strength/D*(a1/(a1+D)) + self.J*drive_strength**3*a1**2*(3*a1**3+11*a1**2*D+15*a1*D**2+9*D**3)/2/D**3/(a1+D)**3/(a1+2*D)/(3*a1+2*D)

    def ix_strength(self, drive_strength):
        D = self.w1 - self.w2
        a1, a2 = self.a1, self.a2
        self.ix = -self.J*drive_strength/(D + a1) + D*a1*self.J*drive_strength**3 / (D + a1)**3 / (2*D + a1) / (2*D + 3*a1)

    def get_transmon_object(self, C):
        def tune_frequency(x, C):
            EC = e_charge**2/(C)/2/planck_h * 1e-9 / 2 / np.pi
            tmon = scq.Transmon(EJmax=x, EC=EC, ng=0, ncut=11)
            return abs(tmon.E01() - self.wc)
        EJ = fsolve(tune_frequency, 10, args=(self.C1))
        tmon_tuned = scq.Transmon(EJmax=EJ, EC=C, ng=0, ncut=11)
        return tmon_tuned
    
    def get_coupled_system(self):
        self.tmon1 = self.get_transmon_object(self.C1)
        self.tmon2 = self.get_transmon_object(self.C2)
        cavity = scq.Oscillator(E_osc=self.wc, truncated_dim=11)
        hilbertspace = HilbertSpace([self.tmon1, self.tmon2, cavity])
        zpf_correction1 = (self.tmon1.EJ/8/self.tmon1.EC)**(1/4)*1/np.sqrt(2)
        zpf_correction2 = (self.tmon2.EJ/8/self.tmon2.EC)**(1/4)*1/np.sqrt(2)
        g1_k = self.g1 / zpf_correction1
        g2_k = self.g2 / zpf_correction2
        hilbertspace.add_interaction(
            g = g1_k,
            op1 = self.tmon1.n_operator,
            op2 = cavity.creation_operator,
            add_hc = True
            )

        hilbertspace.add_interaction(
            g = g2_k,
            op1 = self.tmon2.n_operator,
            op2 = cavity.creation_operator,
            add_hc = True
        )
        hilbertspace.generate_lookup()
        self.hilbertspace = hilbertspace

    def diagonalize(self):
        self.get_coupled_system()
        E_01 = self.hilbertspace.energy_by_bare_index((0,1,0), subtract_ground = True)  # Qubit 1 in ground state, qubit 2 in excited state (|01⟩)
        E_10 = self.hilbertspace.energy_by_bare_index((1,0,0), subtract_ground = True)  # Qubit 1 in excited state, qubit 2 in ground state (|10⟩)
        E_11 = self.hilbertspace.energy_by_bare_index((1,1,0), subtract_ground = True)  # Both qubits in excited state (|11⟩)

        E_02 = self.hilbertspace.energy_by_bare_index((0,2,0), subtract_ground = True)
        E_20 = self.hilbertspace.energy_by_bare_index((2,0,0), subtract_ground = True)

        # Calculate the ZZ interaction
        self.zz_static_numeric = (E_11 - E_10 - E_01)
        self.w1_dressed = E_10
        self.w2_dressed = E_01
        self.a1_dressed = E_20 - 2*E_10
        self.a2_dressed = E_02 - 2*E_01
        self.a1 = self.tmon1.anharmonicity()
        self.a2 = self.tmon2.anharmonicity()
