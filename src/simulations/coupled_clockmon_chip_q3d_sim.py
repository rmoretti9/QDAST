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
from qdast.simulations.coupled_clockmons_q3d_sim import TwoClockmonsQ3DSim

sim_tool = "q3d"
# Simulation parameters
sim_class = TwoClockmonsQ3DSim

sim_parameters = {
    "name": "coupled_clockmons_chip",
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(5e3, 5e3)),
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
        "ansys_tool": sim_tool,
        "path": dir_path,
        "exit_after_run": True,
        "percent_error": 0.1,
        "maximum_passes": 25,
        "minimum_passes": 2,
        "minimum_converged_passes": 3,
        "post_process": PostProcess("produce_cmatrix_table.py"),
    }
)

# Get layout
layout = get_active_or_new_layout()

# Sweep simulations
simulations = [sim_class(layout, **sim_parameters)]
oas = export_simulation_oas(simulations, dir_path)

open_with_klayout_or_default_application(oas)
export_ansys(simulations, **ansys_export_parameters, skip_errors=True)
