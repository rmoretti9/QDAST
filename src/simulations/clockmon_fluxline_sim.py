import logging
import sys
from pathlib import Path

from kqcircuits.qubits.swissmon import Swissmon
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
from qdast.qubits.clockmon import Clockmon
from qdast.dimensioned_chips.two_clockmons_00 import TwoClockmons00

@add_parameters_from(Clockmon)
class ClockmonFluxlineSim(Simulation):
    flux_simulation = Param(pdt.TypeBoolean, "Setup for magnetic flux simulations", True)

    def build(self):
        params = self.get_parameters()
        if not self.flux_simulation:
            params["junction_type"] = "Sim"

        qubit_cell = self.add_element(TwoClockmons00, **params)
        qubit_trans = pya.DTrans(0, False, self.box.center())
        _, refp = self.insert_cell(qubit_cell, qubit_trans, rec_levels=None)

        # Add internal port next to edge
        self.produce_waveguide_to_port(refp["port_flux"], refp["port_flux_corner"], 1, use_internal_ports="at_edge")

        if self.flux_simulation:
            # Setup squid and refinement boxes for magnetic flux simulations
            self.partition_regions = [
                {"name": "squid", "face": "1t1", "vertical_dimensions": 0.0, "region": pya.DBox(995, 840, 1005, 850)},
                {"name": "refine", "face": "1t1", "vertical_dimensions": 5.0, "region": pya.DBox(985, 825, 1015, 855)},
            ]
        else:
            # Add qubit port
            self.ports.append(InternalPort(2, *self.etched_line(refp["port_squid_a"], refp["port_squid_b"])))


# Prepare output directory
dir_path = create_or_empty_tmp_directory(Path(__file__).stem + "_output")

# Simulation parameters
sim_class = ClockmonFluxlineSim  # pylint: disable=invalid-name
sim_parameters = {
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(2000, 2000)),
}

# Get layout
logging.basicConfig(level=logging.WARN, stream=sys.stdout)
layout = get_active_or_new_layout()

# Sweep simulation and solution type
simulations = [
    (
        sim_class(layout, **sim_parameters, name="fluxline_mut_ind", flux_simulation=True),
        AnsysCurrentSolution(
            max_delta_e=0.01,
            frequency=0.1,
            maximum_passes=20,
            integrate_magnetic_flux=True,
            mesh_size={"squid": 2, "vacuumrefine": 4, "substraterefine": 4},
        ),
    ),
    (
        sim_class(layout, **sim_parameters, name="fluxline_decay", flux_simulation=False),
        AnsysHfssSolution(max_delta_s=0.005, frequency=4.5, maximum_passes=20, sweep_enabled=False),
    ),
]

# Export simulation files
export_ansys(simulations, path=dir_path, exit_after_run=True, post_process=PostProcess("calculate_q_from_s.py"))

# Write and open oas file
open_with_klayout_or_default_application(export_simulation_oas(simulations, dir_path))
