from qdast.chips.detection_device_2s_1a import DetectionDevice2s1a
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
import pandas as pd
import os
 

#  feedline_0,737.909,384.738,830.8050000000001,2190.0,50.0,239.0,69.0,300.0
# feedline_1,737.909,384.738,830.8050000000001,2285.0,50.0,239.0,69.0,300.0
# feedline_2,737.909,384.738,830.8050000000001,2120.0,50.0,239.0,69.0,280.0

@add_parameters_from(DetectionDevice2s1a)
class DetectionDevice2s1a00(DetectionDevice2s1a):
    def build(self):
        self.readout_res_lengths = [8223.4, 8371.6, 7995.90]
        self.coupler_res_length = [8720.8, 9023.1]
        self.feedline_meander_lengths_s1 = [2190, 300]
        self.feedline_meander_lengths_s2 = [2285, 300]
        self.feedline_meander_lengths_a = [2120, 280]
        self.n_fingers = [2.4528, 2.4787, 2.4278]
        self.sensing_1_coupler_widths = [205.79, 172.58]
        self.sensing_2_coupler_widths = [230.14, 233.44]
        self.ancilla_coupler_widths = [198.37, 235.46, 198.45]
        self.feedline_capacitor_n_fingers =[3.4432, 3.4432, 3.4432]
        self.name_chip = ""
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