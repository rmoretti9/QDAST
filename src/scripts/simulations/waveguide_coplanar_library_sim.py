import logging
import sys
from pathlib import Path

import numpy as np

from qdast.qubits.clockmon import Clockmon
from qdast.elements.digit_tee_coupler import DigitTee
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.single_element_simulation import (
    get_single_element_sim_class,
)
from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.elmer.elmer_export import export_elmer
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.export.simulation_export import (
    cross_sweep_simulation,
    export_simulation_oas,
)

from qdast.elements.waveguide_segment import WaveguideSegment
sim_tool = "q3d"

# Simulation parameters
sim_class = get_single_element_sim_class(
    WaveguideSegment, #ignore_ports=["port_0", "port_1"]
)  # pylint: disable=invalid-name
sim_parameters = {
    "name": "waveguide",
    "use_internal_ports": True,
    # "use_ports": True,
    "with_squid": False,
    "face_stack": ["1t1"],
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(1000, 1000)),
    # "separate_island_internal_ports": sim_tool != "eigenmode",  # DoublePads specific
    "a": 10,
    "b": 6,
    "waveguide_length": 0,

}

dir_path = create_or_empty_tmp_directory(Path(__file__).stem + f"_{sim_tool}")

# Add eigenmode and Q3D specific settings
# fmt: off
export_parameters_ansys = {
    "ansys_tool": sim_tool,
    "path": dir_path,
    "exit_after_run": True,
    'percent_error': 0.1,
    'maximum_passes': 25,
    'minimum_passes': 2,
    'minimum_converged_passes': 3,
    "post_process": PostProcess("produce_cmatrix_table.py"),
}

# Get layout
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
layout = get_active_or_new_layout()

simulations = []
waveguide_length = []

# Multi face finger (interdigital) capacitor sweeps
for wg_length in np.linspace(10, 500, 101):
    simulations += [
        sim_class(layout, **{**sim_parameters,
        "wg_length": wg_length,
        "name": f"waveguide_coupler_width_{wg_length}",}, #pya.DPath([pya.DPoint(0, 0), pya.DPoint(wg_length, 0)], 0)}
        )
    ]


# Create simulation
oas = export_simulation_oas(simulations, dir_path)

export_ansys(simulations, **export_parameters_ansys)

logging.info(f"Total simulations: {len(simulations)}")
open_with_klayout_or_default_application(oas)
