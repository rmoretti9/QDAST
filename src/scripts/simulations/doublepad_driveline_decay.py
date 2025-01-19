# This paper has an interesting formula to link the Y matrix to the Q-factor, and is consistent with im(Y)/re(Y)
# https://arxiv.org/pdf/2202.06202
# Other ref: https://arxiv.org/pdf/1504.04353
# https://journals.aps.org/prxquantum/pdf/10.1103/PRXQuantum.2.040202

import logging
import sys
from pathlib import Path

from qdast.qubits.clockmon import Clockmon
from kqcircuits.simulations.export.ansys.ansys_solution import AnsysCurrentSolution, AnsysHfssSolution
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
from qdast.dimensioned_chips.single_doublepads_00 import SingleDoublepads00


class TwoClockmonsDrivelineDecay(Simulation):

    def build(self):
        chip = self.add_element(
                    SingleDoublepads00,
                    sim_tool = "q3d",
                    with_squid = False,
                    n = 32          
                )
        self.cell.insert(pya.DCellInstArray(chip.cell_index(), pya.DTrans(0, False, 0, 0)))
        _, refpoints = self.insert_cell(chip)
        # self.ports.append(InternalPort(2, *self.etched_line(refpoints["qb_0_port_island1_signal"], refpoints["qb_0_port_island1_ground"])))
        self.ports.append(InternalPort(3, *self.etched_line(refpoints["qb_0_port_island2_signal"], refpoints["qb_0_port_island2_ground"])))

# Prepare output directory
dir_path = create_or_empty_tmp_directory(Path(__file__).stem + "2_output")

# Simulation parameters
sim_class = TwoClockmonsDrivelineDecay  # pylint: disable=invalid-name
sim_parameters = {
    # "box": pya.DBox(pya.DPoint(0, 200), pya.DPoint(2900, 2200)), #qb1
    "box": pya.DBox(pya.DPoint(1800, 2200), pya.DPoint(4500, 4300)), #qb2
}

# Get layout
logging.basicConfig(level=logging.WARN, stream=sys.stdout)
layout = get_active_or_new_layout()

# Sweep simulation and solution type
simulations = [
    (
        sim_class(layout, **sim_parameters, name="fluxline_decay", flux_simulation=False),
        AnsysHfssSolution(max_delta_s=0.005, frequency=4.5, maximum_passes=20, sweep_enabled=False),
    ),
]

# Export simulation files
export_ansys(simulations, path=dir_path, exit_after_run=False, post_process=PostProcess("calculate_q_from_s.py"))

# Write and open oas file
open_with_klayout_or_default_application(export_simulation_oas(simulations, dir_path))
