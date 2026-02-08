import math

from kqcircuits.elements.element import Element

from kqcircuits.util.parameters import Param, pdt
from kqcircuits.qubits.qubit import Qubit
from kqcircuits.pya_resolver import pya
from kqcircuits.util.refpoints import WaveguideToSimPort, JunctionSimPort
from kqcircuits.util.geometry_helper import arc_points
from qdast.elements.fluxlines.fluxline_tapered import FluxlineTapered


class DoublepadQDAST(Qubit):

    ground_gap = Param(
        pdt.TypeList, "Width, height of the ground gap (µm, µm)", [1000, 900]
    )
    ground_gap_r = Param(pdt.TypeDouble, "Ground gap rounding radius", 50, unit="μm")
    coupler_width = Param(
        pdt.TypeDouble, "Width of the coupler in µm", 150, unit="μm"
    )
    coupler_height = Param(
        pdt.TypeDouble, "Height of the coupler in µm", 20, unit="μm"
    )

    coupler_r = Param(pdt.TypeDouble, "Coupler rounding radius", 5, unit="μm")
    coupler_a = Param(
        pdt.TypeDouble,
        "Width of the coupler waveguide center conductor",
        Element.a,
        unit="μm",
    )
    coupler_offset = Param(
        pdt.TypeDouble, "Distance between coupler and qubit origins", 285, unit="μm"
    )
    drive_position = Param(
        pdt.TypeList, "Coordinate for the drive port (µm, µm)", [-450, 0]
    )
    island_extent = Param(pdt.TypeList, "Islands width and height (µm, µm)", [700, 170])
    island_to_island_distance = Param(
        pdt.TypeDouble, "Island to island distance", 170, unit="μm"
    )

    wire_radius = Param(
        pdt.TypeDouble, "Radius of circular wire region", 75, unit="μm"
    )
    sim_tool = Param(pdt.TypeString, "Simulation tool", "none", choices=["none", "q3d", "eig"])

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
        island1_region = self._build_island()
        island2_region = island1_region.transformed(pya.CplxTrans(1, 0, True, 0, 0))
        
        island_region = island1_region + island2_region
        etch_region = ground_gap_region - island_region
        # Readout coupler
        trans = [(1, 0, False, 0, 0)]
        trans_refp = []
        for t in trans:
            trans_refp.append(t[:-1] + (t[-1]*self.layout.dbu, ))
        coupler_region = pya.Region()
        coupler_gap_region = pya.Region()
        port_directions = [pya.DVector(0, 1)]

        if float(self.coupler_width) != 0:
            coupler_region_add, coupler_gap_region_add = self._build_coupler(0, pya.CplxTrans(*trans[0]),
                                                                                pya.CplxTrans(*trans_refp[0]))
            coupler_region += coupler_region_add
            coupler_gap_region += coupler_gap_region_add
            
            self.add_port(
                    str(0),
                    self.refpoints[str(0)],
                    direction=port_directions[0],
                )
                            
        etch_region += coupler_gap_region
        etch_region -= coupler_region
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

        # Drive port
        self.add_port(
            "drive",
            pya.DPoint(float(self.drive_position[0]), float(self.drive_position[1])),
            direction=pya.DVector(
                float(self.drive_position[0]), float(self.drive_position[1])
            ),
        )


    def _build_coupler(self, coupler_id, trans, trans_refp):
        width = float(self.coupler_width)
        height = float(self.coupler_height)
        offset = float(self.coupler_offset)

        if offset > float(self.ground_gap[1]) / 2 - height - 50:
            has_gap = True
            stem_height = 50
            self.refpoints[str(coupler_id)] = pya.DPoint(0, offset + height + stem_height)*trans_refp
        else:
            has_gap = False
            stem_height = float(self.ground_gap[1]) / 2 - height
            self.refpoints[str(coupler_id)] = pya.DPoint(0, float(self.ground_gap[1]) / 2)*trans_refp

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

    def _build_island(self):
        width = float(self.island_extent[0])
        height = float(self.island_extent[1])
        distance = self.island_to_island_distance

        island_points = [
                    pya.DPoint(width / 2, distance / 2),]
        island_points +=arc_points(height/2, start=-math.pi/2, stop=math.pi/2, origin=pya.DPoint(width / 2, distance / 2 + height/2))

        island_points += [
                    pya.DPoint(width / 2, distance / 2 + height),
                    pya.DPoint(-width / 2, distance / 2 + height),]
        island_points +=arc_points(height/2, start=math.pi/2, stop=3*math.pi/2, origin=pya.DPoint(-width / 2, distance / 2 + height/2))
        island_points += [
                    pya.DPoint(-width / 2, distance / 2),
                ]
        island_points +=arc_points(float(self.wire_radius), start=0, stop=-math.pi, origin=pya.DPoint(0, distance / 2))
        
        island_region = pya.Region(
            pya.DPolygon(island_points).to_itype(self.layout.dbu))
        
        island_region = island_region.round_corners(
            15 / self.layout.dbu,
            15 / self.layout.dbu,
            self.n,
        )
        # Trapezoid

        width_trapezoid_1 = 75
        width_trapezoid_2 = 5
        height_trapezoid = 82
        trapezoid_points = [pya.DPoint(-width_trapezoid_1 / 2, distance / 2),
                            pya.DPoint(width_trapezoid_1 / 2, distance / 2),
                            pya.DPoint(width_trapezoid_2 / 2, distance / 2 - height_trapezoid),
                            pya.DPoint(-width_trapezoid_2 / 2, distance / 2 - height_trapezoid)]
        trapezoid_region = pya.Region(
            pya.DPolygon(trapezoid_points).to_itype(self.layout.dbu))
        island_region += trapezoid_region
        island_region = island_region.round_corners(
            2 / self.layout.dbu,
            2 / self.layout.dbu,
            self.n,
        )
        # island_region = pya.Region(
        #     pya.DPolygon(
        #         [
        #             pya.DPoint(width / 2, distance / 2),
        #             pya.DPoint(width / 2, distance / 2 + height),
        #             pya.DPoint(-width / 2, distance / 2 + height),
        #             pya.DPoint(-width / 2, distance / 2),
        #         ]
        #     ).to_itype(self.layout.dbu)
        # )

        return island_region

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
                pya.DPoint(4/5*island_width / 2, self.island_to_island_distance / 2 + float(self.island_extent[1])),
                pya.DPoint(
                    4/5*island_width / 2, float(self.ground_gap[1])/2 - port_height / 2
                ),
                pya.DPoint(4/5*island_width / 2, -self.island_to_island_distance / 2 - float(self.island_extent[1])),
                pya.DPoint(
                    4/5*island_width / 2,
                    -float(self.ground_gap[1])/2 + port_height / 2,
                ),
            ]
            port_names = [
                "island1_signal",
                "island1_ground",
                "island2_signal",
                "island2_ground",
            ]
            orientations = [90, 90, -90, -90]
            orientations_ports = [-1, -1, 1, 1]
            ports_offset = [1, 0, -1, 0]
            port_corner_x = [0, 0, 0, 0]
            port_corner_y = [1, -1, -1, 1]
            for i, (trans, name) in enumerate(zip(translations, port_names)):
                port_poly = pya.DCplxTrans(1, orientations[i], False, trans) * port_polygon
                port_regions += pya.Region(port_poly.to_itype(self.layout.dbu))
                self.add_port(
                    name,
                    trans
                    + pya.DPoint(orientations_ports[i]* port_width, ports_offset[i]*port_height / 2),
                    pya.DVector(port_corner_x[i], port_corner_y[i]),
                )

        elif self.sim_tool == "eig":
            port_polygon = pya.DBox(0, 0, self.width_tapered, 4)
            # Eigenmode ports
            leads_length = (
                self.lead_height_untapered
                + self.lead_height_tapered
                + self.bent_section_length
            )
            if self.clock_diameter > self.island_to_island_distance:
                translations = [
                    pya.DPoint(-self.width_tapered/2 + self.external_leads_offset, - leads_length + self.clock_diameter/2),
                    pya.DPoint(-self.width_tapered/2 + self.external_leads_offset, + leads_length - 4 - self.clock_diameter/2),
                ]
            else:
                translations = [
                    pya.DPoint(-self.width_tapered/2 + self.external_leads_offset, - leads_length + self.island_to_island_distance/2),
                    pya.DPoint(-self.width_tapered/2 + self.external_leads_offset, + leads_length - 4 - self.island_to_island_distance/2),
                ]  
            port_names = ["island1", "island2"]
            for i, (trans, name) in enumerate(zip(translations, port_names)):
                port_poly = pya.DCplxTrans(1, 0, False, trans) * port_polygon
                port_regions += pya.Region(port_poly.to_itype(self.layout.dbu))
                self.add_port(
                    name,
                    trans + pya.DPoint(self.width_tapered/2, 4 * (i % 2)),
                    pya.DVector(0, 2 * (i % 2) - 1),
                )

        return port_regions

    def produce_fluxline(self):
        if self.fluxline_type == "none":
            return

        cell = self.add_element(FluxlineTapered)

        refpoints_so_far = self.get_refpoints(self.cell)
        if self.external_leads_offset > 0: 
            trans = pya.DCplxTrans(1, 90, False, self.external_leads_offset + 18.33, 13)
        else:
            trans = pya.DCplxTrans(1, -90, False, self.external_leads_offset - 18.33, -13)
        cell_inst, _ = self.insert_cell(cell, trans)
        self.copy_port("flux", cell_inst)


    @classmethod
    def get_sim_ports(cls, simulation):
        return [
            WaveguideToSimPort("port_0", side="top"),
            JunctionSimPort("port_island1_signal", "port_island1_ground"),
            JunctionSimPort("port_island2_signal", "port_island2_ground"),
        ]
