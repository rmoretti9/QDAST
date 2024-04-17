import logging
import sys
import ast

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
from kqcircuits.simulations.simulation import Simulation
from qdast.qubits.clockmon import Clockmon

from qdast.simulations.single_clockmons_sim import SingleClockmonsSim

sim_tool = "eigenmode"
Lj = 1.43580331e-8

mer_dims_default = {
    "gap": 1.0,
    "vacuum": 1.0,
    "metal": 1.0,
    "substrate": 1.0,
}

# Simulation parameters
sim_class = SingleClockmonsSim

sim_parameters = {
    "name": "clockmon_chip",
    "use_internal_ports": True,
    "use_ports": True,
    "qubit_face": ["1t1"],
    "face_stack": ["1t1"],
    "with_squid": False,
    "sim_tool": "eig",
    "junction_inductance": Lj,
    "box": pya.DBox(pya.DPoint(200, 200), pya.DPoint(2350, 2850)),
    "tls_layer_thickness": [
        4.8e-9 * 1e6,
        0.3e-9 * 1e6,
        2.4e-9 * 1e6,
    ],  # MA, MS, and SA,
    "tls_sheet_approximation": True,
    "metal_height": 0.2,
    "tls_layer_material": ["oxideMA", "oxideMS", "oxideSA"],
    "material_dict": {
        **ast.literal_eval(Simulation.material_dict),
        "oxideMA": {"permittivity": 8},
        "oxideMS": {"permittivity": 11.4},
        "oxideSA": {"permittivity": 4},
    },
}
sim_class.junction_inductance = Lj  # Manually adjusting Lj
sim_class.sim_tool = "eig"

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
        "max_delta_f": 0.1,
        "maximum_passes": 25,
        "minimum_passes": 1,
        "minimum_converged_passes": 2,
        "n_modes": 2,
        "frequency": 1,
        # "mesh_size": {
        #     "1t1_substratemer": 30,
        #     "1t1_vacuummer": 30,
        # },
        "post_process": PostProcess(
            "produce_epr_table.py",
            sheet_approximations={
                "MA": {"thickness": 4.8e-9, "eps_r": 8, "background_eps_r": 1.0},
                "MS": {
                    "thickness": 0.3e-9,
                    "eps_r": 11.4,
                    "background_eps_r": 11.45,
                },
                "SA": {"thickness": 2.4e-9, "eps_r": 4, "background_eps_r": 11.45},
            },
            groups=["MA", "MS", "SA", "substrate", "vacuum"],
        ),
    }
)

#####################
# GEOMETRY SETTINGS #
#####################

# Get layout
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
layout = get_active_or_new_layout()
sweep_parameters_list = [
    {
        "name": "clockmon_chip",  # required when specifying sims manually
        "partition_regions": [
            {
                "name": "mer",
                "face": "1t1",
                "metal_edge_dimensions": [
                    mer_dims_default["gap"],
                    mer_dims_default["metal"],
                ],
                "vertical_dimensions": [
                    mer_dims_default["substrate"],
                    mer_dims_default["vacuum"],
                ],
            },
        ],
    }
]

# Sweep simulations
simulations = []
simulations += [
    sim_class(layout, **{**sim_parameters, **sweep_parameters})
    for sweep_parameters in sweep_parameters_list
]
oas = export_simulation_oas(simulations, dir_path)

open_with_klayout_or_default_application(oas)
export_ansys(simulations, **ansys_export_parameters, skip_errors=True)
