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
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
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
        "port_0",
        "port_1",
        "port_2",
        "port_3",
        "port_4",
        "port_5",
    ],
)
sim_parameters = {
    "name": "clockmon",
    "with_squid": False,
    "face_stack": ["1t1"],
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(1500, 1500)),
    "waveguide_length": 0,
}

dir_path = create_or_empty_tmp_directory(Path(__file__).stem + f"_{sim_tool}")

# Add Q3D specific settings
export_parameters_ansys = {
    "ansys_tool": sim_tool,
    "path": dir_path,
    "exit_after_run": True,
    "percent_error": 0.05,
    "maximum_passes": 25,
    "minimum_passes": 2,
    "minimum_converged_passes": 2,
    "post_process": PostProcess("produce_cmatrix_table.py"),
}

# Get layout
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
layout = get_active_or_new_layout()

simulations = []

name = sim_parameters["name"]
simulations = [
    sim_class(
        layout,
        **{
            **sim_parameters,
            "ground_gap": [630, 610],
            "a": 10,
            "b": 6,
            "island_extent": [535, 200],
            "coupler_widths": [0, 0, 0, 0, 0, 0],
            "island_to_island_distance": 50,
            "coupler_offsets": [0, 0, 0, 0, 0, 0],
            "clock_diameter": 50,
            "bending_angle": 0,
            "sim_tool": "q3d",
            "with_squid": False,
            "pad_width": 6,
            "taper_width": 95 / 7,
            "bent_section_length": 8,
            "lead_height_untapered": 4,
            "lead_height_tapered": 8,
            "drive_position": [0, -405],
        },
    )
]

# Create simulation
oas = export_simulation_oas(simulations, dir_path)

export_ansys(simulations, **export_parameters_ansys)

logging.info(f"Total simulations: {len(simulations)}")
open_with_klayout_or_default_application(oas)
