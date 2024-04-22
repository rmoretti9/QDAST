# This code is part of KQCircuits
# Copyright (C) 2022 IQM Finland Oy
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# https://www.gnu.org/licenses/gpl-3.0.html.
#
# The software distribution should follow IQM trademark policy for open-source software
# (meetiqm.com/developers/osstmpolicy). IQM welcomes contributions to the code. Please see our contribution agreements
# for individuals (meetiqm.com/developers/clas/individual) and organizations (meetiqm.com/developers/clas/organization).

import logging
import sys
from pathlib import Path

import numpy as np

from kqcircuits.qubits.double_pads import DoublePads
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.single_element_simulation import get_single_element_sim_class
from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.elmer.elmer_export import export_elmer
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)


sim_tool = "elmer"

# Simulation parameters
sim_class = get_single_element_sim_class(DoublePads)  # pylint: disable=invalid-name
sim_parameters = {
    "name": "double_pads",
    "use_internal_ports": True,
    "use_ports": True,
    "face_stack": ["1t1"],
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(2000, 2000)),
    "separate_island_internal_ports": sim_tool != "eigenmode",  # DoublePads specific
    "tls_layer_thickness": 5e-3 if sim_tool == "eigenmode" else 0.0,  # in µm
    "tls_sheet_approximation": sim_tool == "eigenmode",
    "waveguide_length": 200,
}

dir_path = create_or_empty_tmp_directory(Path(__file__).stem + f"_output_{sim_tool}")
print(dir_path)
export_parameters_elmer = {
    "tool": "capacitance",
    "workflow": {
        "python_executable": "python",
        "n_workers": 4,
        "elmer_n_processes": 1,
        "gmsh_n_threads": 1,
        "elmer_n_threads": 1,
    },
    "mesh_size": {
        "global_max": 50.0,
        "1t1_gap&1t1_signal": [2.0, 4.0],
        "1t1_gap&1t1_ground": [2.0, 4.0],
    },
}

# Get layout
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
layout = get_active_or_new_layout()

name = sim_parameters["name"]
simulations = [
    sim_class(
        layout,
        **{
            **sim_parameters,
            "ground_gap": [900, 900],
            "a": 5,
            "b": 20,
            "coupler_a": 5,
            "coupler_extent": [50, 20],
            "coupler_offset": 100,
            "junction_type": "Manhattan",
            "island2_taper_junction_width": 31.7,
            "junction_total_length": 39.5,
            "name": name,
        },
    )
]

# Create simulation
oas = export_simulation_oas(simulations, dir_path)

export_elmer(simulations, dir_path, **export_parameters_elmer)

logging.info(f"Total simulations: {len(simulations)}")
open_with_klayout_or_default_application(oas)
