"""
PCell implementation for a SingleDoublepads chip (multiple-qubit layouts).
"""

from math import pi

import pandas as pd
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.elements.meander import Meander
from qdast.qubits.clockmon import Clockmon
from kqcircuits.elements.waveguide_coplanar import WaveguideCoplanar

from kqcircuits.elements.waveguide_coplanar_splitter import (
    WaveguideCoplanarSplitter,
    t_cross_parameters,
)
from kqcircuits.pya_resolver import pya
from kqcircuits.util.coupler_lib import cap_params
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
from kqcircuits.junctions.junction import Junction
from kqcircuits.elements.launcher import Launcher
from kqcircuits.elements.waveguide_coplanar_taper import WaveguideCoplanarTaper
from kqcircuits.util.geometry_helper import arc_points


def _get_num_meanders(meander_length, turn_radius, meander_min_width):
    """Return the integer number of meanders required for a given length.

    The formula matches the original implementation: available straight
    length divided by the effective per-meander width (including turns).
    """
    return int(
        (meander_length - turn_radius * (pi - 2))
        / (meander_min_width + turn_radius * (pi - 2))
    )


@add_parameters_from(
    QDASTChip,
    produce_labels=True,
    with_grid=False,
    marker_types="",
    name_mask="",
    name_chip="",
    name_brand="",
    name_copy="",
)
@add_parameters_from(Junction, "junction_type")
@add_parameters_from(Clockmon, "sim_tool", "with_squid")
@add_parameters_from(QDASTChip, label_offset=1500)
class SingleDoublepads(QDASTChip):
    """PCell declaration for a SingleDoublepads chip.

    This class creates launchers, qubits, readout resonators and feedline
    geometry. Default parameters mirror the original source.
    """

    readout_res_lengths = Param(
        pdt.TypeList,
        "Readout resonator lengths (six resonators)",
        [4000, 4100, 4200, 4300, 4400, 4500],
    )
    n_fingers = Param(
        pdt.TypeList,
        "Number of fingers for readout resonator couplers",
        [1.1, 1.2, 1.3, 1.4, 1.5, 1.6],
    )
    feedline_capacitor_n_fingers = Param(
        pdt.TypeDouble,
        "Number of fingers for feedline resonator coupler",
        3.1,
    )

    type_coupler = Param(
        pdt.TypeList,
        "Coupler type for test resonator couplers",
        ["smooth", "smooth", "smooth", "smooth", "smooth", "smooth"],
    )
    coupler_widths = Param(
        pdt.TypeList,
        "Qubit coupler width",
        [150, 150, 150, 150, 150],
    )
    _readout_structure_info = {
        "feedline": [],
        "tees": [],
        "readout_res_lengths": [],
    }
    l_fingers = Param(
        pdt.TypeList,
        "Length of fingers for readout resonator couplers",
        [30, 30, 30, 30, 30, 30],
    )
    with_feedline_resonator = Param(pdt.TypeBoolean, "Enable feedline resonator", True)

    def build(self):
        """Top-level assembly invoked by the PCell framework.

        Creates launcher elements, places qubits and then assembles the
        feedline, readout resonators and charge-drive lines.
        """
        self.launchers = self.produce_launchers(
            "6-ports-10x10", launcher_assignments={}
        )

        launcher_cell = self.add_element(
            Launcher,
            s=380,
            l=200,
            a_launcher=200,
            b_launcher=154,
            launcher_frame_gap=20,
        )

        # Insert launcher instances and keep key reference points used later
        _, self.cl_2 = self.insert_cell(
            launcher_cell, pya.DCplxTrans(1, 90, False, 453 + 200, 10000 - 620 - 200)
        )
        _, self.cl_3 = self.insert_cell(
            launcher_cell,
            pya.DCplxTrans(1, 90, False, 10000 - 453 - 200, 10000 - 615 - 200),
        )
        _, self.refpoints_fl_out = self.insert_cell(
            launcher_cell, pya.DCplxTrans(1, 0, False, 10000 - 620 - 200, 5000)
        )
        _, self.cl_1 = self.insert_cell(
            launcher_cell, pya.DCplxTrans(1, -90, False, 10000 - 453 - 200, +615 + 200)
        )
        _, self.cl_0 = self.insert_cell(
            launcher_cell, pya.DCplxTrans(1, -90, False, 453 + 200 + 7, +615 + 200)
        )
        _, self.refpoints_fl_in = self.insert_cell(
            launcher_cell, pya.DCplxTrans(1, -180, False, 620 + 200, 5000)
        )

        self._produce_qubits()

        self._produce_feedline_resonator()
        self._produce_readout_resonators()
        self._produce_chargelines()

    def _produce_waveguide(
        self, path, term2=0, turn_radius=None, a=None, b=None, object=None
    ):
        """Add a WaveguideCoplanar element following the provided path.

        When ``object == 'feedline'`` the waveguide length is recorded in
        ``_readout_structure_info['feedline']``.
        """
        if a is None:
            a = self.a
        if b is None:
            b = self.b
        if turn_radius is None:
            turn_radius = self.r
        waveguide = self.add_element(
            WaveguideCoplanar,
            path=pya.DPath(path, 1),
            r=turn_radius,
            term2=term2,
            a=a,
            b=b,
        )
        self.insert_cell(waveguide)
        if object == "feedline":
            self._readout_structure_info["feedline"].append(waveguide.length())
        return waveguide

    def _produce_qubit(self, coupler_width, center_x, center_y, rotation, name=None):
        """Instantiate a Clockmon qubit at the requested position.

        Returns the absolute reference points created by the insertion.
        """
        qubit = self.add_element(
            Clockmon,
            ground_gap=[630, 610],
            a=10,
            b=6,
            island_extent=[535, 200],
            coupler_widths=[coupler_width, 0, 0, 0, 0, 0],
            island_to_island_distance=50,
            coupler_offsets=[255, 0, 0, 0, 0, 0],
            clock_diameter=50,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered",
            sim_tool=self.sim_tool,
            with_squid=self.with_squid,
            pad_width=6,
            taper_width=95 / 7,
            bent_section_length=8,
            lead_height_untapered=4,
            lead_height_tapered=8,
            drive_position=[0, -355],
        )
        qubit_trans = pya.DTrans(rotation, False, center_x, center_y)
        _, refpoints_abs = self.insert_cell(qubit, qubit_trans, name, rec_levels=None)
        return refpoints_abs

    def _produce_qubits(self):
        """Create six qubits arranged around the feedline.

        Coordinates and spacing choices follow the original layout.
        """
        qubit_spacing_x = 2000
        qubit_spacing_x_alt = 700
        qubit_spacing_y = 5e3
        qubits_center_x = 5e3
        y_a = 5e3 + 1500
        y_b = 5e3 - 1500

        qb0_refpoints = self._produce_qubit(
            float(self.coupler_widths[0]),
            qubits_center_x - qubit_spacing_x,
            y_b,
            0,
            "qb_0",
        )
        qb1_refpoints = self._produce_qubit(
            float(self.coupler_widths[1]),
            qubits_center_x - qubit_spacing_x + qubit_spacing_x_alt,
            y_a,
            2,
            "qb_1",
        )
        qb2_refpoints = self._produce_qubit(
            float(self.coupler_widths[2]), qubits_center_x, y_b, 0, "qb_2"
        )
        qb3_refpoints = self._produce_qubit(
            float(self.coupler_widths[3]),
            qubits_center_x + qubit_spacing_x_alt,
            y_a,
            2,
            "qb_3",
        )
        qb4_refpoints = self._produce_qubit(
            float(self.coupler_widths[0]),
            qubits_center_x + qubit_spacing_x,
            y_b,
            0,
            "qb_4",
        )
        qb5_refpoints = self._produce_qubit(
            float(self.coupler_widths[1]),
            qubits_center_x + qubit_spacing_x + qubit_spacing_x_alt,
            y_a,
            2,
            "qb_5",
        )

        self._qubit_x_coords = [
            qubits_center_x - qubit_spacing_x,
            qubits_center_x - qubit_spacing_x + qubit_spacing_x_alt,
            qubits_center_x,
            qubits_center_x + qubit_spacing_x_alt,
            qubits_center_x + qubit_spacing_x,
            qubits_center_x + qubit_spacing_x_alt + qubit_spacing_x,
        ]
        self._qubit_refpoints = [
            qb0_refpoints,
            qb1_refpoints,
            qb2_refpoints,
            qb3_refpoints,
            qb4_refpoints,
            qb5_refpoints,
        ]

    def _produce_readout_resonator(
        self, capacitor, capacitor_dtrans, res_idx, qubit_refp
    ):
        """Create a readout resonator connected to a capacitive coupler and a tee.

        A meander section is inserted to match the requested resonator length.
        """
        flip = -1 if res_idx % 2 else 1
        total_length = float(self.readout_res_lengths[res_idx])
        turn_radius = 50
        w = 500
        capacitor_pos = self.get_refpoints(capacitor, capacitor_dtrans)["port_a"]

        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a
            ),
        )
        cross_length = self.b * 3 + 2 * self.a
        cross_trans1 = pya.CplxTrans(
            1,
            -90,
            False,
            pya.DPoint(
                qubit_refp["port_0"].x,
                capacitor_pos.y - flip * (150 + self.a / 2 + self.b),
            ),
        )
        _, tee_refpoints = self.insert_cell(cell_cross, cross_trans1)

        tee_port_1 = "port_right" if res_idx % 2 else "port_left"
        tee_port_2 = "port_left" if res_idx % 2 else "port_right"
        wg_1 = self._produce_waveguide(
            [capacitor_pos, pya.DPoint(capacitor_pos.x, tee_refpoints[tee_port_1].y)]
        )
        wg_2 = self._produce_waveguide(
            [
                tee_refpoints[tee_port_2],
                pya.DPoint(qubit_refp["port_0"].x, qubit_refp["port_0"].y + flip * 300),
                pya.DPoint(qubit_refp["port_0"].x, qubit_refp["port_0"].y),
            ]
        )
        wg_3 = self._produce_waveguide(
            [
                tee_refpoints["port_bottom"],
                pya.DPoint(
                    tee_refpoints["port_bottom"].x - 300, tee_refpoints["port_bottom"].y
                ),
                pya.DPoint(
                    tee_refpoints["port_bottom"].x - 500, tee_refpoints["port_bottom"].y
                ),
                pya.DPoint(
                    tee_refpoints["port_bottom"].x - 500,
                    tee_refpoints["port_bottom"].y - flip * 100,
                ),
            ],
            turn_radius=turn_radius,
        )

        meander_length = (
            total_length - wg_1.length() - cross_length - wg_2.length() - wg_3.length()
        )
        num_meanders = _get_num_meanders(meander_length, turn_radius, w)
        self.insert_cell(
            Meander,
            start_point=pya.DPoint(
                tee_refpoints["port_bottom"].x - 500,
                tee_refpoints["port_bottom"].y - flip * 100,
            ),
            end_point=pya.DPoint(
                tee_refpoints["port_bottom"].x - 500,
                tee_refpoints["port_bottom"].y - flip * 750,
            ),
            length=meander_length,
            meanders=num_meanders,
            r=turn_radius,
        )

    def _produce_feedline_resonator(self):
        """Assemble a feedline resonator with multiple tees and a feedline capacitor.

        The routine places six tees, creates a coupling capacitor and routes
        the waveguides to connect in/out reference points.
        """
        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a
            ),
        )
        tee_refpoints = []
        tee_rotations = [0, 2, 0, 2, 0, 2]
        for i in range(6):
            cross_trans = pya.DTrans(
                tee_rotations[i], False, self._qubit_refpoints[i]["port_0"].x, 5e3
            )
            inst_cross, _ = self.insert_cell(cell_cross, cross_trans)
            tee_refpoints.append(self.get_refpoints(cell_cross, inst_cross.dtrans))
        self._readout_structure_info["tees"] = [self.a, self.a, 2 * self.a]
        self._tee_refpoints = tee_refpoints

        feedline_tee_trans = pya.DTrans(0, False, 8600, 5e3)
        feedline_tee, feedline_tee_refp = self.insert_cell(
            cell_cross, feedline_tee_trans
        )

        if isinstance(self.feedline_capacitor_n_fingers, tuple):
            self.feedline_capacitor_n_fingers = self.feedline_capacitor_n_fingers[0]
        cplr_params = cap_params(
            self.feedline_capacitor_n_fingers,
            float(self.l_fingers[0]),
            "smooth",
            finger_gap=self.b,
        )
        cplr = self.add_element(**cplr_params)
        cplr_refpoints_rel = self.get_refpoints(cplr)
        cplr_dtrans = pya.DTrans(
            0,
            False,
            self.refpoints_fl_in["base"].x + 200,
            self.refpoints_fl_in["base"].y,
        )
        _, cplr_refpoints = self.insert_cell(cplr, cplr_dtrans, rec_levels=None)

        self._produce_waveguide(
            [self.refpoints_fl_in["base"], cplr_refpoints["port_a"]], object="feedline"
        )
        self._produce_waveguide(
            [
                cplr_refpoints["port_b"],
                pya.DPoint(
                    cplr_refpoints["port_b"].x + 200, cplr_refpoints["port_b"].y
                ),
                pya.DPoint(cplr_refpoints["port_b"].x + 200, 5200),
                pya.DPoint(cplr_refpoints["port_b"].x + 200, 5950),
                pya.DPoint(cplr_refpoints["port_b"].x + 400, 5950),
                pya.DPoint(cplr_refpoints["port_b"].x + 400, 5000),
                pya.DPoint(cplr_refpoints["port_b"].x + 600, 5000),
                pya.DPoint(cplr_refpoints["port_b"].x + 600, 5950),
                pya.DPoint(cplr_refpoints["port_b"].x + 800, 5950),
                pya.DPoint(cplr_refpoints["port_b"].x + 800, 5e3),
                tee_refpoints[0]["port_left"],
            ],
            object="feedline",
        )

        # Chain the remaining tees and connect to the out reference points
        self._produce_waveguide(
            [tee_refpoints[0]["port_right"], tee_refpoints[1]["port_right"]],
            object="feedline",
        )
        self._produce_waveguide(
            [tee_refpoints[1]["port_left"], tee_refpoints[2]["port_left"]],
            object="feedline",
        )
        self._produce_waveguide(
            [tee_refpoints[2]["port_right"], tee_refpoints[3]["port_right"]],
            object="feedline",
        )
        self._produce_waveguide(
            [tee_refpoints[3]["port_left"], tee_refpoints[4]["port_left"]],
            object="feedline",
        )
        self._produce_waveguide(
            [tee_refpoints[4]["port_right"], tee_refpoints[5]["port_right"]],
            object="feedline",
        )
        self._produce_waveguide(
            [tee_refpoints[5]["port_left"], feedline_tee_refp["port_left"]],
            object="feedline",
        )
        self._produce_waveguide(
            [feedline_tee_refp["port_right"], self.refpoints_fl_out["base"]],
            object="feedline",
        )
        self._produce_waveguide(
            [
                feedline_tee_refp["port_bottom"],
                pya.DPoint(
                    feedline_tee_refp["port_bottom"].x,
                    feedline_tee_refp["port_bottom"].y - 500,
                ),
                pya.DPoint(
                    feedline_tee_refp["port_bottom"].x + 200,
                    feedline_tee_refp["port_bottom"].y - 500,
                ),
                pya.DPoint(
                    feedline_tee_refp["port_bottom"].x + 200,
                    feedline_tee_refp["port_bottom"].y - 1200,
                ),
            ],
            object="feedline",
        )

    def _produce_readout_resonators(self):
        """Create readout resonators for each coupler/tee pair."""
        for i in range(6):
            cplr_params = cap_params(
                float(self.n_fingers[i]),
                float(self.l_fingers[i]),
                self.type_coupler[i],
                finger_gap=self.b,
            )
            cplr = self.add_element(**cplr_params)
            cplr_refpoints_rel = self.get_refpoints(cplr)
            if i % 2 == 0:
                cplr_pos = (
                    self._tee_refpoints[i]["port_bottom"]
                    - pya.DTrans.R90 * cplr_refpoints_rel["port_b"]
                )
            else:
                cplr_pos = (
                    self._tee_refpoints[i]["port_bottom"]
                    + pya.DTrans.R90 * cplr_refpoints_rel["port_b"]
                )
            cplr_dtrans = pya.DTrans(2 * (i % 2) + 1, False, cplr_pos.x, cplr_pos.y)
            self.insert_cell(cplr, cplr_dtrans)

            self._produce_readout_resonator(
                cplr, cplr_dtrans, i, self._qubit_refpoints[i]
            )

    def _produce_chargelines(self):
        """Create charge-drive lines and tapers connecting launchers to qubit drives.

        The method places two tees, creates tapers for impedance transition
        and routes the necessary waveguide segments for both top and bottom
        rows of qubits.
        """
        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a
            ),
        )
        cross_length = self.b * 3 + 2 * self.a
        cross_trans_0 = pya.CplxTrans(
            1,
            90,
            False,
            pya.DPoint(
                self._qubit_refpoints[0]["port_drive"].x,
                self.cl_0["base"].y + 1000 + self.b + self.a / 2,
            ),
        )
        cross_trans_1 = pya.CplxTrans(
            1,
            -90,
            False,
            pya.DPoint(
                self._qubit_refpoints[5]["port_drive"].x,
                self.cl_2["base"].y - 1000 - self.b - self.a / 2,
            ),
        )
        _, tee_refpoints_0 = self.insert_cell(cell_cross, cross_trans_0)
        _, tee_refpoints_1 = self.insert_cell(cell_cross, cross_trans_1)

        taper_cell, _ = WaveguideCoplanarTaper.create_with_refpoints(
            self.layout,
            self.LIBRARY_NAME,
            a=self.a,
            b=self.b,
            a2=self.a / 3,
            b2=self.b / 3,
            taper_length=80,
        )
        _, taper_ref0 = self.insert_cell(
            taper_cell,
            pya.CplxTrans(
                1,
                90,
                False,
                pya.DPoint(
                    tee_refpoints_0["port_right"].x,
                    tee_refpoints_0["port_right"].y + 1000,
                ),
            ),
        )
        _, taper_ref1 = self.insert_cell(
            taper_cell,
            pya.CplxTrans(
                1,
                -90,
                False,
                pya.DPoint(
                    self._qubit_refpoints[1]["port_drive"].x,
                    tee_refpoints_1["port_right"].y - 1000,
                ),
            ),
        )
        _, taper_ref2 = self.insert_cell(
            taper_cell,
            pya.CplxTrans(
                1,
                90,
                False,
                pya.DPoint(
                    self._qubit_refpoints[2]["port_drive"].x,
                    tee_refpoints_0["port_right"].y + 1000,
                ),
            ),
        )
        _, taper_ref3 = self.insert_cell(
            taper_cell,
            pya.CplxTrans(
                1,
                -90,
                False,
                pya.DPoint(
                    self._qubit_refpoints[3]["port_drive"].x,
                    tee_refpoints_1["port_right"].y - 1000,
                ),
            ),
        )
        _, taper_ref4 = self.insert_cell(
            taper_cell,
            pya.CplxTrans(
                1,
                90,
                False,
                pya.DPoint(
                    self._qubit_refpoints[4]["port_drive"].x,
                    tee_refpoints_0["port_right"].y + 1000,
                ),
            ),
        )
        _, taper_ref5 = self.insert_cell(
            taper_cell,
            pya.CplxTrans(
                1,
                -90,
                False,
                pya.DPoint(
                    self._qubit_refpoints[5]["port_drive"].x,
                    tee_refpoints_1["port_right"].y - 1000,
                ),
            ),
        )

        # Tapered lines connecting tapers to qubit drive ports
        self._produce_waveguide(
            [taper_ref0["port_b"], self._qubit_refpoints[0]["port_drive"]],
            a=self.a / 3,
            b=self.b / 3,
            term2=self.b,
        )
        self._produce_waveguide(
            [taper_ref1["port_b"], self._qubit_refpoints[1]["port_drive"]],
            a=self.a / 3,
            b=self.b / 3,
            term2=self.b,
        )
        self._produce_waveguide(
            [taper_ref2["port_b"], self._qubit_refpoints[2]["port_drive"]],
            a=self.a / 3,
            b=self.b / 3,
            term2=self.b,
        )
        self._produce_waveguide(
            [taper_ref3["port_b"], self._qubit_refpoints[3]["port_drive"]],
            a=self.a / 3,
            b=self.b / 3,
            term2=self.b,
        )
        self._produce_waveguide(
            [taper_ref4["port_b"], self._qubit_refpoints[4]["port_drive"]],
            a=self.a / 3,
            b=self.b / 3,
            term2=self.b,
        )
        self._produce_waveguide(
            [taper_ref5["port_b"], self._qubit_refpoints[5]["port_drive"]],
            a=self.a / 3,
            b=self.b / 3,
            term2=self.b,
        )

        # BOTTOM ROW WIRING
        self._produce_waveguide(
            [
                self.cl_0["base"],
                pya.DPoint(self.cl_0["base"].x, self.cl_0["base"].y + 200),
                pya.DPoint(
                    self._qubit_refpoints[0]["port_drive"].x, self.cl_0["base"].y + 200
                ),
                pya.DPoint(
                    self._qubit_refpoints[0]["port_drive"].x, self.cl_0["base"].y + 300
                ),
            ]
        )
        self._produce_waveguide(
            [
                pya.DPoint(
                    self._qubit_refpoints[0]["port_drive"].x, self.cl_0["base"].y + 300
                ),
                tee_refpoints_0["port_left"],
            ]
        )
        self._produce_waveguide(
            [
                tee_refpoints_0["port_right"],
                pya.DPoint(
                    tee_refpoints_0["port_right"].x,
                    tee_refpoints_0["port_right"].y + 1000,
                ),
            ]
        )

        self._produce_waveguide(
            [
                tee_refpoints_0["port_bottom"],
                pya.DPoint(
                    self._qubit_refpoints[2]["port_drive"].x,
                    tee_refpoints_0["port_bottom"].y,
                ),
                taper_ref2["port_a"],
            ]
        )
        self._produce_waveguide(
            [
                self.cl_1["base"],
                pya.DPoint(self.cl_1["base"].x, self.cl_1["base"].y + 200),
                pya.DPoint(
                    self._qubit_refpoints[4]["port_drive"].x, self.cl_1["base"].y + 200
                ),
                pya.DPoint(
                    self._qubit_refpoints[4]["port_drive"].x, taper_ref4["port_a"].y
                ),
            ]
        )

        # TOP ROW WIRING
        self._produce_waveguide(
            [
                self.cl_2["base"],
                pya.DPoint(self.cl_2["base"].x, self.cl_2["base"].y - 200),
                pya.DPoint(
                    self._qubit_refpoints[1]["port_drive"].x, self.cl_2["base"].y - 200
                ),
                pya.DPoint(
                    self._qubit_refpoints[1]["port_drive"].x, taper_ref1["port_a"].y
                ),
            ]
        )
        self._produce_waveguide(
            [
                self.cl_3["base"],
                pya.DPoint(self.cl_3["base"].x, self.cl_3["base"].y - 200),
                pya.DPoint(
                    self._qubit_refpoints[5]["port_drive"].x, self.cl_3["base"].y - 200
                ),
                pya.DPoint(
                    self._qubit_refpoints[5]["port_drive"].x, self.cl_3["base"].y - 300
                ),
            ]
        )
        self._produce_waveguide(
            [
                tee_refpoints_1["port_right"],
                pya.DPoint(
                    tee_refpoints_1["port_right"].x,
                    tee_refpoints_1["port_right"].y - 1000,
                ),
            ]
        )
        self._produce_waveguide(
            [
                tee_refpoints_1["port_bottom"],
                pya.DPoint(
                    self._qubit_refpoints[3]["port_drive"].x,
                    tee_refpoints_1["port_bottom"].y,
                ),
                pya.DPoint(
                    self._qubit_refpoints[3]["port_drive"].x, taper_ref3["port_a"].y
                ),
            ]
        )
        self._produce_waveguide(
            [
                pya.DPoint(
                    self._qubit_refpoints[5]["port_drive"].x, self.cl_3["base"].y - 300
                ),
                tee_refpoints_1["port_left"],
            ]
        )
