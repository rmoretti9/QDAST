from qdast.chips.detection_device_2s_1a import DetectionDevice2s1a
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os
 
@add_parameters_from(DetectionDevice2s1a)
class DetectionDevice2s1a00(DetectionDevice2s1a):
    def build(self):
        self.readout_res_lengths = [9000, 7000, 8200]
        self.feedline_meander_lengths_s1 = [1980, 300]
        self.feedline_meander_lengths_s2 = [1980, 300]
        self.feedline_meander_lengths_a = [1980, 300]
        self.n_fingers = [1.8884, 1.5, 1.3]
        self.sensing_1_coupler_widths = [100, 100]
        self.sensing_2_coupler_widths = [100, 100]
        self.ancilla_coupler_widths = [150, 80, 90]
        self.feedline_capacitor_n_fingers =[3.4432, 4.3, 4]
        self.name_chip = "V00"
        self.name_copy = ""
        self.with_squid = False
        self.margin = 120
        super().build()
        self.get_readout_structure_info()

    def get_readout_structure_info(self):
        self._readout_structure_info["readout_structure"] = ["feedline_resonator"]

        df = pd.DataFrame.from_dict(self._readout_structure_info, orient = "index")
        current_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_directory, 'detection_device_2s_1a_00.csv')
        df.to_csv(file_path, mode='w')