import logging
import sys
from pathlib import Path

from qdast.qubits.clockmon import Clockmon
from kqcircuits.simulations.export.ansys.ansys_solution import AnsysCurrentSolution, AnsysEigenmodeSolution
from kqcircuits.simulations.port import InternalPort
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.simulation import Simulation
from kqcircuits.util.parameters import add_parameters_from, Param, pdt

from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from qdast.dimensioned_chips.single_doublepads_04 import SingleDoublepads04


class SingleDoublepadsSim(Simulation):

    def build(self):
        chip = self.add_element(
                    SingleDoublepads04,
                    sim_tool = "q3d",
                    with_squid = False,
                    n = 32          
                )
        self.cell.insert(pya.DCellInstArray(chip.cell_index(), pya.DTrans(0, False, 0, 0)))

# Prepare output directory
dir_path = create_or_empty_tmp_directory(Path(__file__).stem + "_output")

# Simulation parameters
sim_class = SingleDoublepadsSim  # pylint: disable=invalid-name
sim_parameters = {
    "box": pya.DBox(pya.DPoint(200, 200), pya.DPoint(7500-200, 7300)),
}
ansys_export_parameters = {
    "path": dir_path,
    "exit_after_run": False,
}
ansys_export_parameters.update(
    {
        "ansys_tool": "eigenmode",
        "max_delta_f": 0.1,  # quite tight
        "maximum_passes": 3,
        "minimum_passes": 1,
        "minimum_converged_passes": 2,
        "n_modes": 1,
        "min_frequency": 3,  # minimum allowed frequency
        # mer_correction_path is not portable! TODO: make relative path
        "post_process": PostProcess(
            "produce_epr_table.py",
        ),
    }
)
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