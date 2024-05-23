from qdast.chips.single_clockmons import SingleClockmons
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd

@add_parameters_from(SingleClockmons)
class SingleClockmons01(SingleClockmons):
    def build(self):
        self.readout_res_lengths = [8220.1, 8101.4, 7986.1, 7873.9]
        self.n_fingers = [2.0865, 2.0687, 2.0513, 2.0344]
        self.coupler_widths = [127.82, 131.24, 134.53, 137.69]
        self.with_feedline_resonator = True
        self._readout_structure_info = {
        "feedline": [],
        "tees": [],
        "readout_res_lengths": [],
    }
        super().build()
        self.get_readout_structure_info()

    def get_readout_structure_info(self):
        if self.with_feedline_resonator:
            self._readout_structure_info["readout_structure"] = ["feedline_resonator"]
        else:
            self._readout_structure_info["readout_structure"] = ["simple"]

        df = pd.DataFrame.from_dict(self._readout_structure_info, orient = "index")
        df.to_csv("C:/Users/labranca/Desktop/work/QDAST/src/modeling/single_clockmons/single_clockmons_readout_structure01.csv", mode='w')
