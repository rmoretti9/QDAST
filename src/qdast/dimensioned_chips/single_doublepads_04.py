from qdast.chips.single_doublepads_75_75 import SingleDoublepads7575
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os


@add_parameters_from(SingleDoublepads7575)
class SingleDoublepads04(SingleDoublepads7575):
    def build(self):
        self.readout_res_lengths = [3990.8, 3887.0, 3782.5, 3680.0]
        self.n_fingers = [1.8398, 1.7783, 1.767, 1.7831]
        self.coupler_widths = [291.76, 276.6, 262.61, 249.69]
        self.with_feedline_resonator = True
        self.feedline_capacitor_n_fingers = 4.25267767
        self.feedline_meander_height = 1350
        self._readout_structure_info = {
            "feedline": [],
            "tees": [],
            "readout_res_lengths": [],
        }
        self.name_chip = "V04"
        self.name_copy = ""
        self.with_squid = False
        self.margin = 120
        super().build()
        self.get_readout_structure_info()

    def get_readout_structure_info(self):
        if self.with_feedline_resonator:
            self._readout_structure_info["readout_structure"] = ["feedline_resonator"]
        else:
            self._readout_structure_info["readout_structure"] = ["simple"]

        df = pd.DataFrame.from_dict(self._readout_structure_info, orient="index")
        current_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_directory, "single_doublepads_04.csv")
        df.to_csv(file_path, mode="w")
