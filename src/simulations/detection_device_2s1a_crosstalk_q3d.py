from qdast import defaults
from kqcircuits import defaults
from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from qdast.chips.detection_device_2s_1a import DetectionDevice2s1a
from kqcircuits.simulations.simulation import Simulation


class DetectionDevice2S1ASim(Simulation):

    def build(self):
        chip = self.add_element(
            DetectionDevice2s1a, sim_tool="q3d", with_squid=False, n=24
        )
        # self.cell.insert(pya.DCellInstArray(chip.cell_index(), pya.DTrans(0, False, 0, 0)))
        _, refpoints = self.insert_cell(chip)


sim_tool = "q3d"
# Simulation parameters
sim_class = DetectionDevice2S1ASim

sim_parameters = {
    "name": "detection_device_2s1a",
    "box": pya.DBox(pya.DPoint(2500, 3000), pya.DPoint(8000, 8000)),
}
dir_path = create_or_empty_tmp_directory(f"{sim_parameters['name']}_{sim_tool}")

#######################
# SIMULATION SETTINGS #
#######################

ansys_export_parameters = {
    "path": dir_path,
    "exit_after_run": False,
}

ansys_export_parameters.update(
    {
        "ansys_tool": sim_tool,
        "path": dir_path,
        "percent_error": 0.1,
        "maximum_passes": 25,
        "minimum_passes": 2,
        "minimum_converged_passes": 3,
        "post_process": PostProcess("produce_cmatrix_table.py"),
    }
)

# Get layout
layout = get_active_or_new_layout()

# Sweep simulations
simulations = [sim_class(layout, **sim_parameters)]
oas = export_simulation_oas(simulations, dir_path)

open_with_klayout_or_default_application(oas)
export_ansys(simulations, **ansys_export_parameters, skip_errors=True)
