# This paper has an interesting formula to link the Y matrix to the Q-factor, and is consistent with im(Y)/re(Y)
# https://arxiv.org/pdf/2202.06202
# Other ref: https://arxiv.org/pdf/1504.04353
# https://journals.aps.org/prxquantum/pdf/10.1103/PRXQuantum.2.040202

import logging
import sys
from pathlib import Path

from qdast.qubits.clockmon import Clockmon
from kqcircuits.simulations.export.ansys.ansys_solution import (
    AnsysCurrentSolution,
)
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.simulation import Simulation

from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from qdast.qubits.clockmon import Clockmon


class ClockmonFluxlineSim(Simulation):

    def build(self):
        chip = self.add_element(
            Clockmon,
            sim_tool="q3d",
            with_squid=False,
            n=32,
            ground_gap=[630, 690],
            a=10,
            b=6,
            island_extent=[540, 240],
            coupler_widths=[100, 0, 0, 100, 0, 0],
            island_to_island_distance=50,
            coupler_offsets=[295, 0, 0, 295, 0, 0],
            clock_diameter=0,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered",
            pad_width=6,
            taper_width=0,
            width_tapered=4,
            width_untapered=4,
            drive_position=[0, 0],
            external_leads_offset=305,
            bent_section_length=6,
            lead_height_tapered=4,
            lead_height_untapered=4,
            fluxline_type="Fluxline Tapered",
        )
        self.cell.insert(
            pya.DCellInstArray(chip.cell_index(), pya.DTrans(0, False, 0, 0))
        )
        _, refpoints = self.insert_cell(chip)
        self.produce_waveguide_to_port(
            refpoints["port_flux"],
            refpoints["port_flux_corner"],
            1,
            use_internal_ports="at_edge",
        )

        self.partition_regions = [
            {
                "name": "squid",
                "face": "1t1",
                "vertical_dimensions": 0.0,
                "region": pya.DBox(297, -8, 312, 8),
            },
            {
                "name": "refine",
                "face": "1t1",
                "vertical_dimensions": 5.0,
                "region": pya.DBox(280, -20, 330, 20),
            },
        ]


# Prepare output directory
dir_path = create_or_empty_tmp_directory(Path(__file__).stem + "_output")

# Simulation parameters
sim_class = ClockmonFluxlineSim  # pylint: disable=invalid-name
sim_parameters = {
    "box": pya.DBox(pya.DPoint(-1000, -1000), pya.DPoint(1000, 1000)),
    "port_island_1": True,
    "port_island_2": False,
}

# Get layout
logging.basicConfig(level=logging.WARN, stream=sys.stdout)
layout = get_active_or_new_layout()

# Sweep simulation and solution type
simulations = [
    (
        sim_class(
            layout, **sim_parameters, name="fluxline_mut_ind", flux_simulation=True
        ),
        AnsysCurrentSolution(
            max_delta_e=0.01,
            frequency=0.1,
            maximum_passes=20,
            integrate_magnetic_flux=True,
            mesh_size={"squid": 2, "vacuumrefine": 4, "substraterefine": 4},
        ),
    ),
]

# Export simulation files
export_ansys(
    simulations,
    path=dir_path,
    exit_after_run=False,
    post_process=PostProcess("calculate_q_from_s.py"),
)

# Write and open oas file
open_with_klayout_or_default_application(export_simulation_oas(simulations, dir_path))
