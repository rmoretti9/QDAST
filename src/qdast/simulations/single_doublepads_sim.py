from kqcircuits.simulations.simulation import Simulation
from kqcircuits.pya_resolver import pya
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
from qdast.dimensioned_chips.single_doublepads_00 import SingleDoublepads00
from qdast.chips.single_clockmons import SingleClockmons

from kqcircuits.simulations.port import InternalPort
from kqcircuits.elements.fluxlines.fluxline import Fluxline
from kqcircuits.junctions.junction import Junction

class SingleDoublepadsSim(Simulation):

    def build(self):
        chip = self.add_element(
                    SingleDoublepads00,
                    sim_tool = "eig",
                    with_squid = False,
                    n = 24          
                )
        self.cell.insert(pya.DCellInstArray(chip.cell_index(), pya.DTrans(0, False, 0, 0)))
        _, refpoints = self.insert_cell(chip)
        self.ports.extend(
            [
                InternalPort(
                    1,
                    *self.etched_line(refpoints["qb_0_port_island1"], refpoints["qb_0_port_island2"]),
                    inductance=self.junction_inductance,
                    capacitance=1e-16,
                    junction=True,
                ),])
            # self.ports