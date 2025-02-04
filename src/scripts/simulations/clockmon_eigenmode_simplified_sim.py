import logging
import sys
import ast

from qdast import defaults
from kqcircuits import defaults
from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.elmer.elmer_export import export_elmer
from kqcircuits.simulations.export.xsection.xsection_export import (
    create_xsections_from_simulations,
    separate_signal_layer_shapes,
)
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
    get_simulation_directory,
)
from kqcircuits.simulations.simulation import Simulation
from kqcircuits.simulations.single_element_simulation import (
    get_single_element_sim_class,
)
from qdast.qubits.clockmon import Clockmon
from kqcircuits.simulations.export.xsection.xsection_export import (
    create_xsections_from_simulations,
    separate_signal_layer_shapes,
    visualise_xsection_cut_on_original_layout,
)
from kqcircuits.simulations.export.xsection.xsection_export import (
    create_xsections_from_simulations,
    separate_signal_layer_shapes,
    visualise_xsection_cut_on_original_layout,
)

# sim_tools = ["elmer", "eigenmode"] # eigenmode
sim_tool = "eigenmode"
Lj = 1.4348926834765362e-08

# Simulation parameters
sim_class = get_single_element_sim_class(
    Clockmon,
    ignore_ports=["port_drive", "port_1", "port_2", "port_3", "port_4", "port_5", "port_island1_signal", "port_island2_signal"],
)
sim_class.junction_inductance = Lj # Manually adjusting Lj
sim_class.sim_tool = "eig"
sim_parameters = {
    "name": f"clockmon",
    # 'use_internal_ports': True,
    # 'use_ports': True,
    "qubit_face": ["1t1"],
    "face_stack": ["1t1"],
    "with_squid": False,
    "junction_inductance": Lj,
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(2000, 2000)),
    "tls_sheet_approximation": True,
    "ground_gap": [630, 610],
    "a": 10,
    "b": 6,
    "coupler_widths": [150, 0, 0, 0, 0, 0],
    "island_extent": [535, 200], 
    "island_to_island_distance": 50,
    "coupler_offsets": [255, 0, 0, 0, 0, 0],
    "clock_diameter": 90,
    "sim_tool": "eig",
    "bending_angle": 0,
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
        "ansys_tool": "eigenmode",
        "max_delta_f": 0.05,  # quite tight
        "maximum_passes": 25,
        "minimum_passes": 1,
        "minimum_converged_passes": 2,
        "n_modes": 1,
        "min_frequency": 1,  # minimum allowed frequency
        "mesh_size": {
            "1t1_substratemer": 20,
            "1t1_vacuummer": 20,
        },
        # mer_correction_path is not portable! TODO: make relative path
        "post_process": PostProcess(
            "produce_epr_table.py",
        ),
    }
)

   
    #####################
    # GEOMETRY SETTINGS #
    #####################

# Get layout
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
layout = get_active_or_new_layout()
    # layout.dbu = 1e-5  # need finer DBU for SA-interface

metal_edge_region_x = 40.0
sweep_parameters_list = [
    {
        "name": "clockmon",  # required when specifying sims manually

    }
]

# Sweep simulations
simulations = []  # [sim_class(layout, **sim_parameters)]
simulations += [
    sim_class(layout, **{**sim_parameters, **sweep_parameters})
    for sweep_parameters in sweep_parameters_list
]

    
export_ansys(simulations, **ansys_export_parameters, skip_errors=True)

