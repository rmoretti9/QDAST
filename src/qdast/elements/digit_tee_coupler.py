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

class DigitTee(SmoothCapacitor):
    waveguide_length = Param(pdt.TypeDouble, "Waveguide_length", 200, unit="μm")

    def build(self):
        smooth_capacitor = self.add_element(SmoothCapacitor)
        _, refp = self.insert_cell(smooth_capacitor)
        
        tee = self.add_element(
            WaveguideCoplanarSplitter,
            lengths = [150, 150, 30]
            )
        tee_trans = pya.DTrans(
            1, False, refp['port_a'].x - 30, refp['port_a'].y
        )
        
        _, refp_tee = self.insert_cell(tee, tee_trans)
        
        wg = self.add_element(
            WaveguideCoplanar,
            path = pya.DPath([pya.DPoint(0, 0), pya.DPoint(self.waveguide_length, 0)], 0)
        ) 
        wg_trans = pya.DTrans(
            0, False, refp['port_b'].x, refp['port_b'].y
        )
        _, refp_wg = self.insert_cell(wg, wg_trans)

        # Coupler port
        self.add_port(
            "tee_signal",
            refp_tee["port_a"],
            direction=pya.DVector(pya.DPoint(0, -1)),
        )
        self.add_port(
            "tee_ground",
            refp_tee["port_a"] + pya.DPoint(0, 12),
            direction=pya.DVector(pya.DPoint(0, 1)),
        )
        self.add_port(
            "wg_signal",
            refp_wg["port_b"],
            direction=pya.DVector(pya.DPoint(1, 0)),
        )
        self.add_port(
            "wg_ground",
            refp_wg["port_b"] + pya.DPoint(12, 0),
            direction=pya.DVector(pya.DPoint(1, 0)),
        )        
        open_to_ground = pya.Region(pya.DBox(-self.a/2 - self.b, 0, self.a/2 + self.b, 12).to_itype(self.layout.dbu))
        open_to_ground_trans = [pya.DTrans(0, False, refp_tee['port_b'].x, refp['port_b'].y + 150),
                                pya.DTrans(0, False, refp_tee['port_a'].x, refp['port_a'].y - 150 -12),
                                pya.DTrans(math.pi/2, False, refp_wg['port_b'].x + 12, refp_wg['port_b'].y)
                                ]

        for t in open_to_ground_trans:
            self.cell.shapes(self.get_layer("base_metal_gap_wo_grid")).insert(open_to_ground, t)

    @classmethod
    def get_sim_ports(cls, simulation):
        return [
            JunctionSimPort("port_tee_signal", "port_tee_ground"),
            JunctionSimPort("port_wg_signal", "port_wg_ground"),
        ]
