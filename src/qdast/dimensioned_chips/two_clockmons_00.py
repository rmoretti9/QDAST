from qdast.chips.two_clockmons import TwoClockmons
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os

@add_parameters_from(TwoClockmons)
class TwoClockmons00(TwoClockmons):
    def build(self):
        self.readout_res_lengths = [3849.03, 3958.36] 
        self.coupler_length = 8279.7
        self.n_fingers = [1.9151, 1.9476]
        self.readout_coupler_widths = [95.018, 93.116]
        self.qubit_coupler_widths = [166.09, 166.09]
        self.with_feedline_resonator = True,
        self.name_chip = "V00"
        self.name_copy = "UniMiB"
        self.with_squid = False
        self.margin = 120
        super().build()
        self.get_readout_structure_info()

    def get_readout_structure_info(self):
        if self.with_feedline_resonator:
            self._readout_structure_info["readout_structure"] = ["feedline_resonator"]
        else:
            self._readout_structure_info["readout_structure"] = ["simple"]

        df = pd.DataFrame.from_dict(self._readout_structure_info, orient = "index")
        current_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_directory, 'two_clockmons_00.csv')
        df.to_csv(file_path, mode='w')