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
from kqcircuits.simulations.export.elmer.elmer_export import export_elmer
from kqcircuits.simulations.export.simulation_export import cross_sweep_simulation, export_simulation_oas
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from kqcircuits.simulations.post_process import PostProcess

sim_tool = "q3d"

# Simulation parameters
sim_class = get_single_element_sim_class(
    Clockmon, ignore_ports=["port_drive", "port_island1", "port_island2", "port_0", "port_1", "port_2"]
)  # pylint: disable=invalid-name

# Get layout
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
layout = get_active_or_new_layout()

sim_parameters = {
    "name": "clockmon",
    # "use_internal_ports": True,
    # "use_ports": True,
    "with_squid": False,
    "face_stack": ["1t1"],
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(1500, 1500)),
    # "separate_island_internal_ports": sim_tool != "eigenmode",  # DoublePads specific
    "waveguide_length": 200,
    "ground_gap": [630, 610],
    "a": 10,
    "b": 6,
    "island_extent": [535, 200],
    "island_to_island_distance": 50,
    "coupler_offsets": [0, 0, 0, 255, 285, 285],
    "coupler_heights": [0, 0, 0, 20, 20, 20],
    "clock_diameter": 0,
    "bending_angle": 0,
    "sim_tool": "q3d",
    "pad_width": 6,
    "taper_width": 0,
    "width_tapered": 4,
    "width_untapered": 6,
    "drive_position": [0, 0],
    "external_leads_offset": 305,
    "bent_section_length": 6,
    "lead_height_tapered": 4,
    "lead_height_untapered": 4,
}


dir_path = create_or_empty_tmp_directory(Path(__file__).stem + f"_{sim_tool}")

# Add eigenmode and Q3D specific settings
# fmt: off
export_parameters_ansys = {
    "ansys_tool": sim_tool,
    "path": dir_path,
    "exit_after_run": True,
    'percent_error': 0.2,
    'maximum_passes': 25,
    'minimum_passes': 2,
    'minimum_converged_passes': 2,
    "post_process": PostProcess("produce_cmatrix_table.py"),
}

simulations = []
coupler_widths_3 = np.linspace(30, 300, 6)
coupler_widths_4 = np.linspace(30, 200, 6)
coupler_widths_5 = np.linspace(30, 200, 6)


coupler_width_set = []
for cp_w_3 in coupler_widths_3:
    for cp_w_4 in coupler_widths_4:
        for cp_w_5 in coupler_widths_5:
            coupler_width_set.append([0, 0, 0, cp_w_3, cp_w_4, cp_w_5])
simulations +=cross_sweep_simulation(
        layout,
        sim_class,
        sim_parameters,
        {"coupler_widths": coupler_width_set}

    )

# Create simulation
oas = export_simulation_oas(simulations, dir_path)

export_ansys(simulations, **export_parameters_ansys)

open_with_klayout_or_default_application(oas)
