from qdast.chips.single_clockmons import SingleClockmons
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd

@add_parameters_from(SingleClockmons)
class SingleClockmons01(SingleClockmons):
    def build(self):
        self.readout_res_lengths = [8294.5, 8172.6 , 8028.1, 7888.3]
        self.n_fingers = [1.9843, 2.0027, 2.1644, 2.2957]
        self.coupler_widths = [112.3, 108.81, 105.43, 102.16]
        self.with_feedline_resonator = True,
        self.feedline_capacitor_n_fingers = 4.25267767,
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
