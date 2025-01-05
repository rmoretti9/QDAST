from qdast.chips.single_doublepads import SingleDoublepads
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os

@add_parameters_from(SingleDoublepads)
class SingleDoublepads00(SingleDoublepads):
    def build(self):
        self.readout_res_lengths = [4200, 4200 , 4100, 4200, 3900, 3800]
        self.n_fingers = [1.9843, 2.0027, 2.1644, 2.2957, 2, 2]
        self.coupler_widths = [112.3, 108.81, 105.43, 102.16, 100, 100]
        self.with_feedline_resonator = True,
        self.feedline_capacitor_n_fingers = 4.25267767,
        self._readout_structure_info = {
            "feedline": [],
            "tees": [],
            "readout_res_lengths": [],
        }
        self.name_chip = "V00"
        self.name_copy = "UniMiB"
        self.with_squid = True
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
        file_path = os.path.join(current_directory, 'single_doublepads_00.csv')
        df.to_csv(file_path, mode='w')