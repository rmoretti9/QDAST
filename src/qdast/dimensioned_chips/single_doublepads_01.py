from qdast.chips.single_doublepads import SingleDoublepads
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os

@add_parameters_from(SingleDoublepads)
class SingleDoublepads01(SingleDoublepads):
    def build(self):
        self.readout_res_lengths = [3841.5, 3774.3 , 3552.6, 3620.7, 3691.4, 3508.5]
        self.n_fingers = [2.4374, 2.4126, 2.3428, 2.3656, 2.3887, 2.328]
        self.coupler_widths = [358.15, 325.3, 306.22, 317.95, 330.52, 298.75]
        self.with_feedline_resonator = True,
        self.feedline_capacitor_n_fingers = 4.25267767,
        self._readout_structure_info = {
            "feedline": [],
            "tees": [],
            "readout_res_lengths": [],
        }
        self.name_chip = "V01"
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

        df = pd.DataFrame.from_dict(self._readout_structure_info, orient = "index")
        current_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_directory, 'single_doublepads_01.csv')
        df.to_csv(file_path, mode='w')