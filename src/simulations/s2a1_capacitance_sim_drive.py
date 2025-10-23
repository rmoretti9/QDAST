# This paper has an interesting formula to link the Y matrix to the Q-factor, and is consistent with im(Y)/re(Y)
# https://arxiv.org/pdf/2202.06202
# Other ref: https://arxiv.org/pdf/1504.04353
# https://journals.aps.org/prxquantum/pdf/10.1103/PRXQuantum.2.040202

import logging
import sys
from pathlib import Path

from kqcircuits.simulations.export.ansys.ansys_solution import (
    AnsysQ3dSolution,
)

from kqcircuits.simulations.simulation import Simulation


from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from qdast.dimensioned_chips.detection_device_2s_1a_00 import DetectionDevice2s1a00


class TwoClockmonsDrivelineDecay(Simulation):

    def build(self):
        chip = self.add_element(
            DetectionDevice2s1a00, sim_tool="q3d", with_squid=False, n=32
        )
        self.cell.insert(
            pya.DCellInstArray(chip.cell_index(), pya.DTrans(0, False, 0, 0))
        )


# Prepare output directory
dir_path = create_or_empty_tmp_directory(Path(__file__).stem + "_output")

# Simulation parameters
sim_class = TwoClockmonsDrivelineDecay  # pylint: disable=invalid-name
sim_parameters = {
    # "box": pya.DBox(pya.DPoint(2000, 3000), pya.DPoint(4800, 5500)), #S1
    # "box": pya.DBox(pya.DPoint(5300, 3000), pya.DPoint(8000, 5500)), #qb1
    "box": pya.DBox(pya.DPoint(3500, 5500), pya.DPoint(6500, 9000)),  # A
}

# Get layout
logging.basicConfig(level=logging.WARN, stream=sys.stdout)
layout = get_active_or_new_layout()

# Sweep simulation and solution type
simulations = [
    (
        sim_class(layout, **sim_parameters),
        AnsysQ3dSolution(
            name="_q3d",
            percent_error=0.02,
            maximum_passes=25,
            minimum_converged_passes=2,
            minimum_passes=2,
        ),
    ),
]

# Export simulation files
export_ansys(simulations, path=dir_path, exit_after_run=False)

# Write and open oas file
open_with_klayout_or_default_application(export_simulation_oas(simulations, dir_path))
