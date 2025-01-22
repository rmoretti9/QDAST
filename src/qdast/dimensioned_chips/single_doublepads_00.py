from qdast.chips.single_doublepads import SingleDoublepads
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os
 
@add_parameters_from(SingleDoublepads)
class SingleDoublepads00(SingleDoublepads):
    def build(self):
        self.readout_res_lengths = [4074.6, 3990.4 , 3757.3, 3831.9, 3909.6, 3709.1]
        self.n_fingers = [1.7994, 1.7762, 1.7127, 1.733, 1.7541, 1.6995]
        self.coupler_widths = [134.7, 128.92, 113.33, 118.27, 123.46, 110.18]
        self.with_feedline_resonator = True,
        self.feedline_capacitor_n_fingers = 4.25267767,
        self._readout_structure_info = {
            "feedline": [],
            "tees": [],
            "readout_res_lengths": [],
        }
        self.name_chip = "V00"
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
        file_path = os.path.join(current_directory, 'single_doublepads_00.csv')
        df.to_csv(file_path, mode='w')