from qdast.chips.single_clockmons import SingleClockmons
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from

@add_parameters_from(SingleClockmons)
class SingleClockmons00(SingleClockmons):
    def build(self):
        self.readout_res_lengths = [4000, 7000, 7000, 7000]
        self.n_fingers = [2.0865, 2.0687, 2.0513, 2.0344]
        super().build()