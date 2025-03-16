from qdast.chips.detection_device_2s_1a import DetectionDevice2s1a
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os
 
@add_parameters_from(DetectionDevice2s1a)
class DetectionDevice2s1a00(DetectionDevice2s1a):
    def build(self):
        self.readout_res_lengths = [7779.67, 7975.48, 7537.83]
        self.coupler_res_length = [8273.20 , 8540.80]
        self.feedline_meander_lengths_s1 = [1980, 300]
        self.feedline_meander_lengths_s2 = [2075, 300]
        self.feedline_meander_lengths_a = [1895, 280]
        self.n_fingers = [2.1464, 2.1804, 2.1132]
        self.sensing_1_coupler_widths = [155.68, 153.38]
        self.sensing_2_coupler_widths = [173.90, 193.11]
        self.ancilla_coupler_widths = [149.3, 154.43, 157.05]
        self.feedline_capacitor_n_fingers =[3.4432, 3.4432, 3.4432]
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