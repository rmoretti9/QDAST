from qdast.chips.single_clockmons import SingleClockmons
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os


@add_parameters_from(SingleClockmons)
class SingleClockmons00(SingleClockmons):
    def build(self):
        self.readout_res_lengths = [8295.0, 8178.9, 8066.0, 7956.1]
        self.n_fingers = [1.9811, 1.9643, 1.9476, 1.9311]
        self.coupler_widths = [112.3, 108.81, 105.43, 102.16]
        self.with_feedline_resonator = False
        self._readout_structure_info = {
            "feedline": [],
            "tees": [],
            "readout_res_lengths": [],
        }
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

        df = pd.DataFrame.from_dict(self._readout_structure_info, orient="index")
        current_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_directory, "single_clockmons_00.csv")
        df.to_csv(file_path, mode="w")
