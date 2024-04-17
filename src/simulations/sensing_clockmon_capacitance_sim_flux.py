import logging
import sys
from pathlib import Path

from qdast.qubits.clockmon import Clockmon
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
)
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.simulations.simulation import Simulation

sim_tool = "q3d"


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
            coupler_widths=[0, 0, 0, 0, 0, 0],
            island_to_island_distance=50,
            coupler_offsets=[0, 0, 0, 0, 0, 0],
            clock_diameter=0,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered",
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


# Simulation parameters
sim_class = ClockmonFluxlineSim

dir_path = create_or_empty_tmp_directory(Path(__file__).stem + f"_{sim_tool}")

# Add eigenmode and Q3D specific settings
# fmt: off
export_parameters_ansys = {
    "ansys_tool": sim_tool,
    "path": dir_path,
    "exit_after_run": False,
    'percent_error': 0.03,
    'maximum_passes': 25,
    'minimum_passes': 2,
    'minimum_converged_passes': 2,
    "post_process": PostProcess("produce_cmatrix_table.py"),
}

# Get layout
logging.basicConfig(level=logging.INFO, stream=sys.stdout)
layout = get_active_or_new_layout()

simulations = []

simulations = [
    sim_class(
        layout,
        box= pya.DBox(pya.DPoint(-1000,-1000), pya.DPoint(1000, 1000)),
    )
]

# Create simulation
oas = export_simulation_oas(simulations, dir_path)

export_ansys(simulations, **export_parameters_ansys)

logging.info(f"Total simulations: {len(simulations)}")
open_with_klayout_or_default_application(oas)
