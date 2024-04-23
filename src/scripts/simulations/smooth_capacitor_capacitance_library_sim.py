import logging
import sys
from pathlib import Path
import numpy as np

from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.simulation_export import (
    cross_sweep_simulation,
    export_simulation_oas,
)

from kqcircuits.elements.smooth_capacitor import SmoothCapacitor
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.single_element_simulation import (
    get_single_element_sim_class,
)
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)

# Prepare output directory
dir_path = create_or_empty_tmp_directory(Path(__file__).stem + "_output")

sim_class = get_single_element_sim_class(
    SmoothCapacitor
)  # pylint: disable=invalid-name

# Simulation parameters, using multiface interdigital as starting point
sim_parameters = {
    "name": "finger_capacitor",
    "use_internal_ports": True,
    "use_ports": True,
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(1000, 1000)),
    "a2": -1,
    "b2": -1,
    "finger_width": 10,
    "finger_gap": 6,
    "ground_padding": 10,
    "port_size": 200,
    "face_stack": ["1t1"],
}
# Parameters that differ from sim_parameters for gap type
gap_parameters = {
    "finger_number": 4,
    "finger_length": 0,
    "finger_gap": 0,
    "finger_width": 10,
}
export_parameters = {
    "path": dir_path,
    "ansys_tool": "q3d",
    "post_process": PostProcess("produce_cmatrix_table.py"),
    "exit_after_run": True,
    "percent_error": 0.2,
    "minimum_converged_passes": 2,
    "maximum_passes": 20,
}

# Sweep ranges
finger_control = np.linspace(1, 8, 2)

# Get layout
logging.basicConfig(level=logging.WARN, stream=sys.stdout)
layout = get_active_or_new_layout()

# Cross sweep number of fingers and finger length
simulations = []

# Multi face finger (interdigital) capacitor sweeps
simulations += cross_sweep_simulation(
    layout,
    sim_class,
    sim_parameters,
    {"finger_control": finger_control},
)

# # Multi face gap capacitor sweeps
# simulations += cross_sweep_simulation(
#     layout,
#     sim_class,
#     {
#         **sim_parameters,
#         "name": sim_parameters["name"] + "_gap",
#         **gap_parameters,
#     },
# )

# Export Ansys files
export_ansys(simulations, **export_parameters)

# Write and open oas file
open_with_klayout_or_default_application(export_simulation_oas(simulations, dir_path))
