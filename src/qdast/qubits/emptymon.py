import math

from kqcircuits.elements.element import Element
from kqcircuits.util.parameters import Param, pdt
from kqcircuits.qubits.qubit import Qubit
from kqcircuits.pya_resolver import pya
from kqcircuits.util.geometry_helper import arc_points

class Emptymon(Qubit):
    """Just gap and a tee coupler
    """

    ground_gap = Param(
        pdt.TypeList, "Width, height of the ground gap (µm, µm)", [700, 700]
    )
    ground_gap_r = Param(pdt.TypeDouble, "Ground gap rounding radius", 50, unit="μm")
    coupler_width = Param(
        pdt.TypeDouble, "Width of the couplers in µm", 150, unit="μm"
    )
    coupler_height = Param(
        pdt.TypeDouble, "Height of the couplers in µm", 20, unit="μm"
    )

    coupler_r = Param(pdt.TypeDouble, "Coupler rounding radius", 5, unit="μm")
    coupler_a = Param(
        pdt.TypeDouble,
        "Width of the coupler waveguide center conductor",
        Element.a,
        unit="μm",
    )
    coupler_offset = Param(
        pdt.TypeDouble, "Distance between coupler and qubit origin", 200, unit="μm"
    )

    drive_position = Param(
        pdt.TypeList, "Coordinate for the drive port (µm, µm)", [-450, 0]
    )

    def build(self):

        # Qubit base
        ground_gap_points = [
            pya.DPoint(float(self.ground_gap[0]) / 2, float(self.ground_gap[1]) / 2),
            pya.DPoint(float(self.ground_gap[0]) / 2, -float(self.ground_gap[1]) / 2),
            pya.DPoint(-float(self.ground_gap[0]) / 2, -float(self.ground_gap[1]) / 2),
            pya.DPoint(-float(self.ground_gap[0]) / 2, float(self.ground_gap[1]) / 2),
        ]
        ground_gap_polygon = pya.DPolygon(ground_gap_points)
        ground_gap_region = pya.Region(ground_gap_polygon.to_itype(self.layout.dbu))
        ground_gap_region.round_corners(
            self.ground_gap_r / self.layout.dbu,
            self.ground_gap_r / self.layout.dbu,
            self.n,
        )

        # Readout coupler
        coupler_region = pya.Region()
        coupler_gap_region = pya.Region()

        if float(self.coupler_width) != 0:
            coupler_region_add, coupler_gap_region_add = self._build_coupler(pya.CplxTrans(1, 0, False, pya.DPoint(0,0)),
                                                                                pya.CplxTrans(1, 0, False, pya.DPoint(0,0)))
            coupler_region += coupler_region_add
            coupler_gap_region += coupler_gap_region_add
            
            self.add_port(
                    "0",
                    self.refpoints["0"],
                    direction=pya.DVector(0, 1),
                )
                            
        etch_region = ground_gap_region + coupler_gap_region - coupler_region

        # Inserting regions
        self.cell.shapes(self.get_layer("base_metal_gap_wo_grid")).insert(etch_region)

        # Protection
        protection_polygon = pya.DPolygon(
            [
                p
                + pya.DVector(
                    math.copysign(self.margin, p.x), math.copysign(self.margin, p.y)
                )
                for p in ground_gap_points
            ]
        )
        protection_region = pya.Region(protection_polygon.to_itype(self.layout.dbu))
        protection_region.round_corners(
            (self.ground_gap_r + self.margin) / self.layout.dbu,
            (self.ground_gap_r + self.margin) / self.layout.dbu,
            self.n,
        )
        self.add_protection(protection_region)

        # Drive port
        self.add_port(
            "drive",
            pya.DPoint(float(self.drive_position[0]), float(self.drive_position[1])),
            direction=pya.DVector(
                float(self.drive_position[0]), float(self.drive_position[1])
            ),
        )


    def _build_coupler(self, trans, trans_refp):
        width = float(self.coupler_width)
        height = float(self.coupler_height)
        offset = float(self.coupler_offset)

        if offset > float(self.ground_gap[1]) / 2 - height - 50:
            has_gap = True
            stem_height = 50
            self.refpoints["0"] = pya.DPoint(0, offset + height + stem_height)*trans_refp
        else:
            has_gap = False
            stem_height = float(self.ground_gap[1]) / 2 - height
            self.refpoints["0"] = pya.DPoint(0, float(self.ground_gap[1]) / 2)*trans_refp

        coupler_points = [
            pya.DPoint(-width / 2, offset + height),
            pya.DPoint(-width / 2, offset),
            pya.DPoint(width / 2, offset),
            pya.DPoint(width / 2, offset + height),
        ]
        
        if float(self.coupler_width) > 2* self.coupler_r:
            coupler_points +=arc_points(
                self.a/2, start=3*math.pi/2, stop=math.pi, origin=pya.DPoint(self.a, offset + height + self.a/2)
            )
            coupler_points += arc_points(
                self.a/2, start=0, stop=-math.pi/2, origin=pya.DPoint(-self.a, offset + height + self.a/2)
            )
        waveguide_points = [
            pya.DPoint(-self.a / 2, offset + stem_height + height),
            pya.DPoint(-self.a / 2, offset + height),
            pya.DPoint(self.a / 2, offset + height),
            pya.DPoint(self.a / 2, offset + stem_height + height),
        ]
        coupler_region = pya.Region(
            pya.DPolygon(coupler_points).to_itype(self.layout.dbu)
        ).round_corners(self.coupler_r / self.layout.dbu, self.coupler_r / self.layout.dbu, self.n)
        waveguide_region = pya.Region(
            pya.DPolygon(waveguide_points).to_itype(self.layout.dbu)
        )
        total_coupler_region = coupler_region + waveguide_region

        if has_gap:
            total_coupler_gap_region = total_coupler_region.sized(
                self.b / self.layout.dbu, self.b / self.layout.dbu, self.n
            )
            etch_top_waveguide_region = pya.Region(
                pya.DBox(
                    -self.a / 2 - self.b,
                    offset + height + stem_height,
                    self.a / 2 + self.b,
                    offset + height + stem_height + self.b,
                ).to_itype(self.layout.dbu)
            )
            total_coupler_gap_region -= etch_top_waveguide_region
            protection_region = total_coupler_gap_region.sized(
                self.margin / self.layout.dbu
            )
            self.add_protection(protection_region.transformed(trans))

        else:
            total_coupler_gap_region = pya.Region()
        total_coupler_region = total_coupler_region.transformed(trans)
        total_coupler_gap_region = total_coupler_gap_region.transformed(trans)
        return total_coupler_region, total_coupler_gap_region

    @classmethod
    def get_sim_ports(cls, simulation):
        return [
        ]
