from qdast.chips.single_doublepads_75_75 import SingleDoublepads7575
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os


@add_parameters_from(SingleDoublepads7575)
class SingleDoublepads07(SingleDoublepads7575):
    def build(self):
        self.readout_res_lengths = [8262.0, 8088.8, 7922.6, 7769.8]
        self.n_fingers = [2.3479, 2.324, 2.2997, 2.2759]
        self.coupler_widths = [160.25, 158.19, 156.18, 144.66]
        self.x_offset = [0, 0, -850, -850]
        self.qubit_type = "doublepad_qdast"
        self.with_feedline_resonator = True
        self.alternate_drivelines = True
        self.feedline_capacitor_n_fingers = 4.25267767
        self.feedline_meander_height = 1350
        self._readout_structure_info = {
            "feedline": [],
            "tees": [],
            "readout_res_lengths": [],
        }
        self.name_chip = "V07"
        self.name_copy = ""
        self.with_squid = False
        self.margin = 120
        self.resonator_type = "half"
        self.tail_variant = "2"
        super().build()
        self.get_readout_structure_info()

    def get_readout_structure_info(self):
        if self.with_feedline_resonator:
            self._readout_structure_info["readout_structure"] = ["feedline_resonator"]
        else:
            self._readout_structure_info["readout_structure"] = ["simple"]

        df = pd.DataFrame.from_dict(self._readout_structure_info, orient="index")
        current_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_directory, "single_doublepads_07.csv")
        df.to_csv(file_path, mode="w")
