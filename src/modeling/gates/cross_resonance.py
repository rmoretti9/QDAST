"""Various calculations to estimate cross-resonance and spurious terms, 
assuming two transmon qubits coupled through a cavity."""
import numpy as np
import scqubits as scq
from scqubits import HilbertSpace
from scipy.optimize import fsolve
import qutip as qt

planck_h = 6.62607015e-34
e_charge = 1.60217663e-19
phi0 =  2.067833848e-15
sx1 = qt.tensor(qt.sigmax(), qt.qeye(2))
sy1 = qt.tensor(qt.sigmay(), qt.qeye(2))
sz1 = qt.tensor(qt.sigmaz(), qt.qeye(2))
sp1 = qt.tensor(qt.sigmap(), qt.qeye(2))
sm1 = qt.tensor(qt.sigmam(), qt.qeye(2))

sx2 = qt.tensor(qt.qeye(2), qt.sigmax())
sy2 = qt.tensor(qt.qeye(2), qt.sigmay())
sz2 = qt.tensor(qt.qeye(2), qt.sigmaz())
sp2 = qt.tensor(qt.qeye(2), qt.sigmap())
sm2 = qt.tensor(qt.qeye(2), qt.sigmam())

class CRModel():
    """Cross-resonance model."""
    def __init__(self, w1, w2, wc, a1, a2, g1, g2, C1, C2, use_dressed_freq, mode_a, truncated_dim, crosstalk = 0):
        self.w1 = w1
        self.w2 = w2
        self.wc = wc
        self.a1_forced = a1
        self.a2_forced = a2
        self.g1 = g1
        self.g2 = g2
        self.C1 = C1
        self.C2 = C2
        self.truncated_dim = truncated_dim

        if not(use_dressed_freq == False and mode_a == "input"):
            self.diagonalize()
        if use_dressed_freq:
            self.w1_use, self.w2_use = self.w1_dressed, self.w2_dressed
        else:
            self.w1_use, self.w2_use = w1, w2
        if mode_a == "numeric_bare":
            self.a1_use, self.a2_use = self.a1_bare, self.a1_bare
        elif mode_a == "numeric_dressed":
            self.a1_use, self.a2_use = self.a1_dressed, self.a2_dressed
        elif mode_a == "input":
            self.a1_use, self.a2_use = a1, a2
        self.compute_J(crosstalk)

    def compute_J(self, crosstalk):
        w1, w2 = self.w1_use, self.w2_use
        numerator = self.g2 * self.g1 * (w1 + w2 - 2 * self.wc)
        denominator = 2 * (w1 - self.wc) * (w2 - self.wc)
        J = numerator / denominator
        self.J = J + crosstalk
    
    def zz_strength(self, drive_strength):
        D = self.w1_use - self.w2_use
        a1, a2 = self.a1_use, self.a2_use
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
        return self.zz_static, self.zz_driven, self.zz
    
    def zx_strength(self, drive_strength):
        D = self.w1_use - self.w2_use
        a1 = self.a1_use
        self.zx = -self.J*drive_strength/D*(a1/(a1+D)) + self.J*drive_strength**3*a1**2*(3*a1**3+11*a1**2*D+15*a1*D**2+9*D**3)/2/D**3/(a1+D)**3/(a1+2*D)/(3*a1+2*D)
        return self.zx

    def ix_strength(self, drive_strength):
        D = self.w1_use - self.w2_use
        a1 = self.a1_use
        self.ix = -self.J*drive_strength/(D + a1) + D*a1*self.J*drive_strength**3 / (D + a1)**3 / (2*D + a1) / (2*D + 3*a1)

    def get_transmon_object(self, w, C):
        def tune_frequency(x, C):
            EC = e_charge**2/(C)/2/planck_h * 1e-9
            tmon = scq.Transmon(EJ=x, EC=EC, ng=0, ncut=11, truncated_dim=self.truncated_dim)
            return abs(tmon.E01() - w)
        EJ = fsolve(tune_frequency, 10, args=(self.C1))[0]
        EC = e_charge**2/(C)/2/planck_h * 1e-9
        tmon_tuned = scq.Transmon(EJ=EJ, EC=EC, ng=0, ncut=11, truncated_dim=self.truncated_dim)
        return tmon_tuned
    
    def get_coupled_system(self):
        self.tmon1 = self.get_transmon_object(self.w1, self.C1)
        self.tmon2 = self.get_transmon_object(self.w2, self.C2)
        self.cavity = scq.Oscillator(E_osc=self.wc, truncated_dim=self.truncated_dim)
        hilbertspace = HilbertSpace([self.tmon1, self.tmon2, self.cavity])
        zpf_correction1 = -1j*(self.tmon1.EJ/8/self.tmon1.EC)**(1/4)*1/np.sqrt(2)
        zpf_correction2 = -1j*(self.tmon2.EJ/8/self.tmon2.EC)**(1/4)*1/np.sqrt(2)
        g1_k = self.g1 / zpf_correction1
        g2_k = self.g2 / zpf_correction2
        hilbertspace.add_interaction(
            g = g1_k,
            op1 = self.tmon1.n_operator,
            op2 = self.cavity.creation_operator,
            add_hc = True
            )

        hilbertspace.add_interaction(
            g = g2_k,
            op1 = self.tmon2.n_operator,
            op2 = self.cavity.creation_operator,
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
        self.zz_static_numeric = (E_11 - E_10 - E_01) / 2
        self.w1_dressed = E_10
        self.w2_dressed = E_01
        self.a1_dressed = E_20 - 2*E_10
        self.a2_dressed = E_02 - 2*E_01
        self.a1_bare = self.tmon1.anharmonicity()
        self.a2_bare = self.tmon2.anharmonicity()

    def time_evolution(self, init_state, drive_strength, tlist, m2, echo=False, tsteps = 101):
        op_list = []
        op_list.extend((sx1, sy1, sz1))
        op_list.extend((sx2, sy2, sz2))
        self.op_list = op_list
        D = self.w1_use - self.w2_use
        Hp = 2*np.pi*(D/2*sz1 + self.J*(sp1*sm2 + sm1*sp2) + drive_strength*sx1 + m2*drive_strength*sx2)
        Hm = 2*np.pi*(D/2*sz1 + self.J*(sp1*sm2 + sm1*sp2) - drive_strength*sx1 - m2*drive_strength*sx2)

        if echo:
            num_operators = len(op_list)
            output = np.zeros((num_operators, len(tlist)))
            for t, time in enumerate(tlist):
                if time == 0:
                    for i in range(num_operators):
                        output[i, :] = qt.expect(op_list[i], init_state)
                    # output(qt.expect(e_op_list, psi0)[0])
                else:
                    time_subarray = np.linspace(0, time, tsteps)
                    output1 = qt.mesolve(Hp, init_state, time_subarray[0:int(len(time_subarray)/2)], [], [])
                    psiEcho = qt.tensor(qt.gates.rx(np.pi), qt.identity(2)) * output1.states[-1]
                    output2 = qt.mesolve(Hm, psiEcho, time_subarray[int(len(time_subarray)/2):], [], op_list)
                    # print(output2.expect[4][-1])
                    for i in range(num_operators):
                        output[i, t] = output2.expect[i][-1]
            return output
        else:
            output = qt.mesolve(Hp, init_state, tlist, [], op_list)
            return output.expect

    def create_report(self):
        print("Frequencies: ", self.w1_use, self.w2_use, "[GHz / 2pi]")
        print("Anharmonicities: ", self.a1_bare, self.a2_bare, "[GHz / 2pi]")
        print("Exchange rate J: ", self.J, "[GHz / 2pi]")