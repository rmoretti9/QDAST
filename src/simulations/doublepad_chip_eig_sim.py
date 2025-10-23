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

import logging
import sys
from qdast.simulations.single_doublepads_sim import SingleDoublepadsSim

Lj = 10.252561e-9
# Simulation parameters
sim_class = SingleDoublepadsSim

sim_parameters = {
    "name": "doublepad_chip",
    "use_internal_ports": True,
    "use_ports": True,
    "qubit_face": ["1t1"],
    "face_stack": ["1t1"],
    "with_squid": False,
    "sim_tool": "eig",
    "junction_inductance": Lj,
    "box": pya.DBox(pya.DPoint(1000, 2200), pya.DPoint(4500, 4200)),
    # "tls_sheet_approximation": True,
}
sim_class.junction_inductance = Lj  # Manually adjusting Lj

dir_path = create_or_empty_tmp_directory(
    f"{sim_parameters['name']}_{sim_parameters['sim_tool']}"
)

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
        "maximum_passes": 3,
        "minimum_passes": 1,
        "minimum_converged_passes": 2,
        "n_modes": 1,
        "frequency": 3,
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
parameters_list = [
    {
        "name": "doublepad_chip",  # required when specifying sims manually
    }
]

# Sweep simulations
simulations = []  # [sim_class(layout, **sim_parameters)]
simulations += [
    sim_class(layout, **{**sim_parameters, **sweep_parameters})
    for sweep_parameters in parameters_list
]
oas = export_simulation_oas(simulations, dir_path)

open_with_klayout_or_default_application(oas)
export_ansys(simulations, **ansys_export_parameters, skip_errors=True)
