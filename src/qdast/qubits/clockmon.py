import math

from kqcircuits.elements.element import Element
from kqcircuits.junctions.squid import Squid
from kqcircuits.junctions.manhattan import Manhattan
from qdast.junctions.manhattan_single_junction_centered import (
    ManhattanSingleJunctionCentered,
)
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
from kqcircuits.qubits.qubit import Qubit
from kqcircuits.pya_resolver import pya
from kqcircuits.util.refpoints import WaveguideToSimPort, JunctionSimPort
from kqcircuits.util.geometry_helper import arc_points


@add_parameters_from(Squid, junction_type="Manhattan Single Junction Centered")
# @add_parameters_from(Manhattan)
@add_parameters_from(ManhattanSingleJunctionCentered)
class Clockmon(Qubit):
    """A two-island qubit, consisting of two rounded rectangles shunted by a junction, with one capacitive coupler.

    Contains a coupler on the north edge and two separate qubit islands in the center
    joined by a junction or SQUID loaded from another library.
    Refpoint for a readout line at the opening to the coupler and a modifiable refpoint for
    a driveline.
    """

    ground_gap = Param(
        pdt.TypeList, "Width, height of the ground gap (µm, µm)", [700, 700]
    )
    ground_gap_r = Param(pdt.TypeDouble, "Ground gap rounding radius", 50, unit="μm")
    coupler_extent = Param(
        pdt.TypeList, "Width, height of the coupler (µm, µm)", [150, 20]
    )
    coupler_r = Param(pdt.TypeDouble, "Coupler rounding radius", 10, unit="μm")
    coupler_a = Param(
        pdt.TypeDouble,
        "Width of the coupler waveguide center conductor",
        Element.a,
        unit="μm",
    )
    coupler_offset = Param(
        pdt.TypeDouble, "Distance between coupler ad qubit origin", 200, unit="μm"
    )
    squid_offset = Param(
        pdt.TypeDouble, "Offset between SQUID center and qubit center", 0, unit="μm"
    )
    drive_position = Param(
        pdt.TypeList, "Coordinate for the drive port (µm, µm)", [-450, 0]
    )
    island_extent = Param(pdt.TypeList, "Islands width and height (µm, µm)", [150, 100])
    island_to_island_distance = Param(
        pdt.TypeDouble, "Island to island distance", 50, unit="μm"
    )
    clock_diameter = Param(pdt.TypeDouble, "Junction space allocation", 70, unit="μm")
    with_squid = Param(pdt.TypeBoolean, "Boolean whether to include the squid", True)
    bending_angle = Param(pdt.TypeDouble, "Leads bending angle", 45, unit="degrees")
    sim_tool = Param(pdt.TypeString, "Simulation tool", "q3d", choices=["q3d", "eig"])

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

        # Islands
        island1_region = self._build_island() + self._build_leads()
        island2_region = island1_region.transformed(pya.CplxTrans(1, 180, False, 0, 0))
        island_region = island1_region + island2_region
        island_region.round_corners(
            20 / self.layout.dbu,
            20 / self.layout.dbu,
            self.n,
        )

        # Readout coupler
        coupler_region, coupler_gap_region = self._build_coupler()

        etch_region = (
            ground_gap_region - island_region + coupler_gap_region - coupler_region
        )

        # SQUID
        squid_transf = pya.DCplxTrans(
            1, -self.bending_angle, False, pya.DVector(0, self.squid_offset)
        )
        if self.with_squid:
            self.produce_squid(
                squid_transf,
                include_base_metal_gap=False,
                include_base_metal_addition=False,
            )
        else:
            port_region = self._make_ports()
            etch_region -= port_region

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

        # Coupler port
        self.add_port(
            "0",
            self.refpoints["cplr"],
            direction=pya.DVector(pya.DPoint(0, float(self.ground_gap[1]))),
        )

        # Drive port
        self.add_port(
            "drive",
            pya.DPoint(float(self.drive_position[0]), float(self.drive_position[1])),
            direction=pya.DVector(
                float(self.drive_position[0]), float(self.drive_position[1])
            ),
        )

    def _build_coupler(self):
        width = float(self.coupler_extent[0])
        height = float(self.coupler_extent[1])
        offset = self.coupler_offset

        if offset > float(self.ground_gap[1]) / 2 - height - 50:
            has_gap = True
            stem_height = 50
            self.refpoints["cplr"] = pya.DPoint(0, offset + height + stem_height)
        else:
            has_gap = False
            stem_height = float(self.ground_gap[1]) / 2 - height
            self.refpoints["cplr"] = pya.DPoint(0, float(self.ground_gap[1]) / 2)

        coupler_points = [
            pya.DPoint(-width / 2, offset + height),
            pya.DPoint(-width / 2, offset),
            pya.DPoint(width / 2, offset),
            pya.DPoint(width / 2, offset + height),
        ]
        waveguide_points = [
            pya.DPoint(-self.a / 2, offset + stem_height + height),
            pya.DPoint(-self.a / 2, offset + height),
            pya.DPoint(self.a / 2, offset + height),
            pya.DPoint(self.a / 2, offset + stem_height + height),
        ]
        coupler_region = pya.Region(
            pya.DPolygon(coupler_points).to_itype(self.layout.dbu)
        ).round_corners(5 / self.layout.dbu, 5 / self.layout.dbu, self.n)
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
            self.add_protection(protection_region)

        else:
            total_coupler_gap_region = pya.Region()

        return total_coupler_region, total_coupler_gap_region

    def _build_island(self):
        width = float(self.island_extent[0])
        height = float(self.island_extent[1])
        distance = self.island_to_island_distance
        island_region = pya.Region(
            pya.DPolygon(
                [
                    pya.DPoint(width / 2, distance / 2),
                    pya.DPoint(width / 2, distance / 2 + height),
                    pya.DPoint(-width / 2, distance / 2 + height),
                    pya.DPoint(-width / 2, distance / 2),
                ]
            ).to_itype(self.layout.dbu)
        )

        if self.clock_diameter > self.island_to_island_distance:
            circular_unecth = pya.Region(
                pya.DPolygon(
                    arc_points(
                        self.clock_diameter / 2,
                        start=math.pi - 0.2,
                        stop=0.2,
                        origin=pya.DPoint(0, 0),
                    )
                ).to_itype(self.layout.dbu)
            )
            island_region -= circular_unecth

        return island_region

    def _build_leads(self):
        self._width_untapered = 8
        self._width_tapered = 4
        self._height_untapered = 10
        self._height_tapered = 8
        self._length_bent_section = 6
        y_offset = (
            self.island_to_island_distance / 2
            if self.island_to_island_distance > self.clock_diameter
            else self.clock_diameter / 2
        )
        bending_angle = self.bending_angle
        lead_points = [
            pya.DPoint(self._width_untapered / 2, y_offset),
            pya.DPoint(-self._width_untapered / 2, y_offset),
            pya.DPoint(-self._width_untapered / 2, y_offset - self._height_untapered),
            pya.DPoint(
                -self._width_tapered / 2,
                y_offset - self._height_untapered - self._height_tapered,
            ),
            pya.DPoint(
                self._width_tapered / 2,
                y_offset - self._height_untapered - self._height_tapered,
            ),
            pya.DPoint(self._width_untapered / 2, y_offset - self._height_untapered),
        ]
        straight_lead_region = pya.Region(
            pya.DPolygon(lead_points).to_itype(self.layout.dbu)
        )
        bent_lead_polygon = pya.DPolygon(
            [
                pya.DPoint(-self._width_tapered / 2, 0),
                pya.DPoint(-self._width_tapered / 2, -self._length_bent_section),
                pya.DPoint(self._width_tapered / 2, -self._length_bent_section),
                pya.DPoint(self._width_tapered / 2, 0),
            ]
        )
        angle_patch_polygon = pya.DPolygon(
            [
                pya.DPoint(-self._width_tapered / 2, 0),
                pya.DPoint(
                    -self._width_tapered / 2 * math.cos(math.radians(bending_angle)),
                    -self._width_tapered / 2 * math.sin(math.radians(bending_angle)),
                ),
                pya.DPoint(0, 0),
            ]
        )
        angle_patch_polygon = (
            pya.DCplxTrans(
                1,
                0,
                False,
                pya.DPoint(0, y_offset - self._height_untapered - self._height_tapered),
            )
            * angle_patch_polygon
        )
        angle_patch_region = pya.Region(angle_patch_polygon.to_itype(self.layout.dbu))

        bent_lead_polygon = (
            pya.DCplxTrans(
                1,
                bending_angle,
                False,
                pya.DPoint(0, y_offset - self._height_untapered - self._height_tapered),
            )
            * bent_lead_polygon
        )
        bent_lead_region = pya.Region(bent_lead_polygon.to_itype(self.layout.dbu))

        lead_region = straight_lead_region + bent_lead_region + angle_patch_region
        return lead_region

    def _make_ports(self):
        port_width = 3
        port_height = 6
        island_width = float(self.island_extent[0])
        gap_width = float(self.ground_gap[0])
        island_midpoint_height = (
            float(self.island_extent[1]) / 2 + self.island_to_island_distance / 2
        )

        port_regions = pya.Region()
        if self.sim_tool == "q3d":
            port_polygon = pya.DBox(0, 0, port_width, port_height)

            # Capacitance matrix ports
            translations = [
                pya.DPoint(island_width / 2, island_midpoint_height - port_height / 2),
                pya.DPoint(
                    gap_width / 2 - port_width, island_midpoint_height - port_height / 2
                ),
                pya.DPoint(island_width / 2, -island_midpoint_height - port_height / 2),
                pya.DPoint(
                    gap_width / 2 - port_width,
                    -island_midpoint_height - port_height / 2,
                ),
            ]
            port_names = [
                "island1_signal",
                "island1_ground",
                "island2_signal",
                "island2_ground",
            ]
            for i, (trans, name) in enumerate(zip(translations, port_names)):
                port_poly = pya.DCplxTrans(1, 0, False, trans) * port_polygon
                port_regions += pya.Region(port_poly.to_itype(self.layout.dbu))
                self.add_port(
                    name,
                    trans
                    + pya.DPoint((-1 * (i % 2) + 1) * port_width, port_height / 2),
                    pya.DVector(-2 * (i % 2) + 1, 0),
                )

        elif self.sim_tool == "eig":
            port_polygon = pya.DBox(0, 0, 4, 2)
            # Eigenmode ports
            leads_length = (
                self._height_untapered
                + self._height_tapered
                + self._length_bent_section
            )
            translations = [
                pya.DPoint(-2, self.clock_diameter / 2 - leads_length),
                pya.DPoint(-2, -(self.clock_diameter / 2 - leads_length + 2)),
            ]
            port_names = ["island1", "island2"]
            for i, (trans, name) in enumerate(zip(translations, port_names)):
                port_poly = pya.DCplxTrans(1, 0, False, trans) * port_polygon
                port_regions += pya.Region(port_poly.to_itype(self.layout.dbu))
                self.add_port(
                    name,
                    trans + pya.DPoint(2, 2 * (i % 2)),
                    pya.DVector(0, 2 * (i % 2) - 1),
                )

        return port_regions

    @classmethod
    def get_sim_ports(cls, simulation):
        return [
            WaveguideToSimPort("port_0", side="top"),
            JunctionSimPort("port_island1_signal", "port_island1_ground"),
            JunctionSimPort("port_island2_signal", "port_island2_ground"),
            JunctionSimPort("port_island1", "port_island2"),
        ]
