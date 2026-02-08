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
    "waveguide_length": 0,
    "sim_tool": sim_tool,
    "use_internal_ports": True
}

dir_path = create_or_empty_tmp_directory(Path(__file__).stem + f"_{sim_tool}")

# Add Q3D specific settings
export_parameters_ansys = {
    "ansys_tool": sim_tool,
    "path": dir_path,
    "exit_after_run": False,
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
            sim_tool: sim_tool,
        },
    )
]

# Create simulation
oas = export_simulation_oas(simulations, dir_path)

export_ansys(simulations, **export_parameters_ansys)

logging.info(f"Total simulations: {len(simulations)}")
open_with_klayout_or_default_application(oas)
