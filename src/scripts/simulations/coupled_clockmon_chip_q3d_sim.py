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
    'percent_error': 0.1,
    'maximum_passes': 25,
    'minimum_passes': 2,
    'minimum_converged_passes': 3,
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

# C1_0 _1 0 1.06053733201e-13
# C1_2 _1 _2 3.66861355889e-14
# C1_3 _1 _3 2.18187498008e-17
# C1_4 _1 _4 1.39279875929e-17
# C2_0 _2 0 1.08514488394e-13
# C2_3 _2 _3 3.29784355009e-17
# C2_4 _2 _4 2.09781223949e-17
# C3_0 _3 0 1.07015490075e-13
# C3_4 _3 _4 3.67196762469e-14
# C4_0 _4 0 1.08496195805e-13