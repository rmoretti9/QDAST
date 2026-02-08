import logging
import sys
from pathlib import Path

import numpy as np

from qdast.qubits.doublepad_qdast import DoublepadQDAST
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.single_element_simulation import (
    get_single_element_sim_class,
)
from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.simulation_export import export_simulation_oas, cross_sweep_simulation
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from kqcircuits.simulations.post_process import PostProcess

sim_tool = "q3d"

# Simulation parameters
sim_class = get_single_element_sim_class(
    DoublepadQDAST,
    ignore_ports=[
        "port_drive",

    ],
)
sim_parameters = {
    "name": "doublepad_qdast",
    "with_squid": False,
    "face_stack": ["1t1"],
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(2500, 2500)),
    "waveguide_length": 55,
    "sim_tool": sim_tool,
    "use_internal_ports": True,
    "ground_gap": [1000, 900],
    "ground_gap_r": 50,
    "coupler_offset" : 285,
    "island_extent": [700, 170],
    "island_to_island_distance": 170,
    "wire_radius": 75,
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


simulations = []
coupler_widths = np.linspace(20, 300, 20)
simulations += cross_sweep_simulation(
    layout,
    sim_class,
    sim_parameters,
    {
        "coupler_width": [
            coupler_width for coupler_width in coupler_widths
        ]
    },
)


# Create simulation
oas = export_simulation_oas(simulations, dir_path)

export_ansys(simulations, **export_parameters_ansys)

open_with_klayout_or_default_application(oas)