import logging
import sys
from pathlib import Path

import numpy as np

from qdast.qubits.clockmon import Clockmon
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.single_element_simulation import (
    get_single_element_sim_class,
)
from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.simulation_export import (
    cross_sweep_simulation,
    export_simulation_oas,
)
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from kqcircuits.simulations.post_process import PostProcess

sim_tool = "q3d"

# Simulation parameters
sim_class = get_single_element_sim_class(
    Clockmon,
    ignore_ports=[
        "port_drive",
        "port_island1",
        "port_island2",
        "port_1",
        "port_2",
        "port_4",
        "port_5",
    ],
)  # pylint: disable=invalid-name

# Get layout
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
layout = get_active_or_new_layout()

sim_parameters = {
    "name": "clockmon",
    "with_squid": False,
    "face_stack": ["1t1"],
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(1500, 1500)),
    "waveguide_length": 200,
    "ground_gap": [630, 610],
    "a": 10,
    "b": 6,
    "island_extent": [535, 200],
    "island_to_island_distance": 50,
    "coupler_offsets": [255, 0, 0, 255, 0, 0],
    "coupler_heights": [20, 0, 0, 20, 0, 0],
    "clock_diameter": 95,
    "bending_angle": 0,
    "sim_tool": "q3d",
    "pad_width": 6,
    "taper_width": 95 / 7,
}

dir_path = create_or_empty_tmp_directory(Path(__file__).stem + f"_{sim_tool}")

# Add Q3D specific settings
export_parameters_ansys = {
    "ansys_tool": sim_tool,
    "path": dir_path,
    "exit_after_run": True,
    "percent_error": 0.2,
    "maximum_passes": 25,
    "minimum_passes": 2,
    "minimum_converged_passes": 2,
    "post_process": PostProcess("produce_cmatrix_table.py"),
}

simulations = []
coupler_widths_0 = np.linspace(20, 200, 21)
coupler_widths_3 = np.linspace(20, 200, 21)

coupler_width_set = []
for cp_w_0 in coupler_widths_0:
    for cp_w_3 in coupler_widths_3:
        coupler_width_set.append([cp_w_0, 0, 0, cp_w_3, 0, 0])
simulations += cross_sweep_simulation(
    layout, sim_class, sim_parameters, {"coupler_widths": coupler_width_set}
)

# Create simulation
oas = export_simulation_oas(simulations, dir_path)

export_ansys(simulations, **export_parameters_ansys)

open_with_klayout_or_default_application(oas)
