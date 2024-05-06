import math

from kqcircuits.elements.element import Element
from kqcircuits.junctions.squid import Squid
from kqcircuits.junctions.manhattan import Manhattan
from qdast.junctions.manhattan_single_junction_centered import (
    ManhattanSingleJunctionCentered,
)
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
from kqcircuits.elements.smooth_capacitor import SmoothCapacitor
from kqcircuits.pya_resolver import pya
from kqcircuits.util.refpoints import WaveguideToSimPort, JunctionSimPort
from kqcircuits.util.geometry_helper import arc_points
from kqcircuits.elements.waveguide_coplanar_splitter import (
    WaveguideCoplanarSplitter,
    t_cross_parameters,
)
from kqcircuits.elements.waveguide_coplanar import WaveguideCoplanar

@add_parameters_from(SmoothCapacitor)
@add_parameters_from(WaveguideCoplanarSplitter)
@add_parameters_from(WaveguideCoplanar)

class WaveguideSegment(SmoothCapacitor):
    wg_length = Param(pdt.TypeDouble, "Waveguide length", 200, unit="μm")

    def build(self):
        
        wg = self.add_element(
            WaveguideCoplanar,
            path = pya.DPath([pya.DPoint(0, 0), pya.DPoint(self.wg_length, 0)], 0)
        ) 
        wg_trans = pya.DTrans(
            0, False, - self.wg_length/2, 0
        )
        _, refp_wg = self.insert_cell(wg, wg_trans)

        self.add_port(
            "0",
            refp_wg["port_a"],
            direction=pya.DVector(pya.DPoint(-1, 0)),
        )
        self.add_port(
            "1",
            refp_wg["port_b"],
            direction=pya.DVector(pya.DPoint(1, 0)),
        )
    @classmethod
    def get_sim_ports(cls, simulation):
        return [
            WaveguideToSimPort("port_0", side="top"),
            WaveguideToSimPort("port_1", side="top"),
        ]
