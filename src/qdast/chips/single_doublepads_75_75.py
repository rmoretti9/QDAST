"""
PCell for a 75x75 SingleDoublepads chip with readout resonators and
feedline geometry.
"""
from math import pi

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
from kqcircuits.elements.waveguide_coplanar_taper import WaveguideCoplanarTaper


def _get_num_meanders(meander_length, turn_radius, meander_min_width):
    """Return the integer number of meanders required for a given length.

    The formula mirrors the original implementation: available straight
    length divided by the per-meander width (including turns).
    """
    return int(
        (meander_length - turn_radius * (pi - 2))
        / (meander_min_width + turn_radius * (pi - 2))
    )


# Add parameters inherited from QDASTChip and other building blocks.
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
class SingleDoublepads7575(QDASTChip):
    """PCell declaration for a 75x75 SingleDoublepads chip.

    Most attributes and methods are preserved verbatim from the original
    implementation; docstrings were added for clarity.
    """

    readout_res_lengths = Param(
        pdt.TypeList,
        "Readout resonator lengths (four resonators)",
        [4000, 4100, 4200, 4300],
    )
    n_fingers = Param(
        pdt.TypeList,
        "Number of fingers for readout resonator couplers",
        [1.1, 1.2, 1.3, 1.4],
    )
    feedline_capacitor_n_fingers = Param(
        pdt.TypeDouble,
        "Number of fingers for feedline resonator coupler",
        3.1,
    )
    feedline_meander_height = Param(
        pdt.TypeDouble,
        "Height of the meandered part of the feedline",
        2100,
    )

    type_coupler = Param(
        pdt.TypeList,
        "Coupler type for test resonator couplers",
        ["smooth", "smooth", "smooth", "smooth"],
    )
    coupler_widths = Param(
        pdt.TypeList,
        "Qubit coupler width",
        [150, 150, 150, 150],
    )
    _readout_structure_info = {
        "feedline": [],
        "tees": [],
        "readout_res_lengths": [],
    }
    l_fingers = Param(
        pdt.TypeList,
        "Length of fingers for readout resonator couplers",
        [30, 30, 30, 30],
    )
    x_offset = Param(
        pdt.TypeList,
        "Readout coupler x-offset",
        [0, 0, 0, 0],
    )
    resonator_type = Param(
        pdt.TypeString, "Readout resonator type", "quarter", choices=["half", "quarter"]
    )
    with_feedline_resonator = Param(
        pdt.TypeBoolean, "Enable feedline resonator", False
    )
    alternate_drivelines = Param(pdt.TypeBoolean, "Alternating drivelines", False)
    tail_variant = Param(pdt.TypeString, "Tail type", "1", choices=["1", "2"])

    def build(self):
        """High-level build sequence invoked by the PCell framework.

        The method creates launchers, qubits, feedline(s) and readout
        resonators according to chosen parameters.
        """
        self.launchers = self.produce_launchers(
            "4-ports-75x75", launcher_assignments={1: "N", 2: "E", 3: "S", 4: "W"}
        )

        self._produce_qubits()

        if not self.with_feedline_resonator:
            self._produce_feedline(200)
        else:
            self._produce_feedline_resonator()
        self._produce_readout_resonators()
        if not self.alternate_drivelines:
            self._produce_chargelines()
        else:
            self._produce_chargelines_v2()

    def _produce_waveguide(self, path, term2=0, turn_radius=None, a=None, b=None, object=None):
        """Helper to add a WaveguideCoplanar element from a sequence of points.

        Returns the created waveguide element and records feedline lengths
        when the optional ``object`` argument equals ``"feedline"``.
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
        """Create a Clockmon qubit and insert it into the layout.

        Returns the absolute reference points produced by the insertion.
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
            with_fluxline=False,
        )
        qubit_trans = pya.DTrans(rotation, False, center_x, center_y)
        _, refpoints_abs = self.insert_cell(qubit, qubit_trans, name, rec_levels=None)
        return refpoints_abs

    def _produce_qubits(self):
        """Create four qubits positioned around the feedline.

        Coordinates and spacing values follow the original layout.
        """
        qubit_spacing_x = 2200
        qubit_spacing_x_alt = 600
        qubit_spacing_y = 3.75e3
        qubits_center_x = 3.75e3

        y_a = 3.75e3 + 1650
        y_b = 3.75e3 - 1650

        qb0_refpoints = self._produce_qubit(
            float(self.coupler_widths[0]), qubits_center_x - qubit_spacing_x + 1200, y_b, 0, "qb_0"
        )
        qb1_refpoints = self._produce_qubit(
            float(self.coupler_widths[1]), qubits_center_x - qubit_spacing_x + qubit_spacing_x_alt + 1200, y_a, 2, "qb_1"
        )
        qb2_refpoints = self._produce_qubit(
            float(self.coupler_widths[2]), qubits_center_x + 1200, y_b, 0, "qb_2"
        )
        qb3_refpoints = self._produce_qubit(
            float(self.coupler_widths[3]), qubits_center_x + qubit_spacing_x_alt + 1200, y_a, 2, "qb_3"
        )

        self._qubit_x_coords = [
            qubits_center_x - qubit_spacing_x,
            qubits_center_x - qubit_spacing_x + qubit_spacing_x_alt,
            qubits_center_x,
            qubits_center_x + qubit_spacing_x_alt,
        ]
        self._qubit_refpoints = [
            qb0_refpoints,
            qb1_refpoints,
            qb2_refpoints,
            qb3_refpoints,
        ]

    def _produce_quarter_readout_resonator(self, capacitor, capacitor_dtrans, res_idx, qubit_refp):
        """Assemble a quarter-wave readout resonator connected via a tee.

        The method creates a splitter, connects waveguides to the capacitor
        and the qubit, and inserts a meander to reach the target length.
        """
        flip = -1 if res_idx % 2 else 1
        total_length = float(self.readout_res_lengths[res_idx])
        turn_radius = 50
        w = 500
        capacitor_pos = self.get_refpoints(capacitor, capacitor_dtrans)["port_a"]

        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a,
            ),
        )
        cross_length = self.b * 3 + 2 * self.a
        cross_trans1 = pya.CplxTrans(
            1,
            -90,
            False,
            pya.DPoint(
                qubit_refp["port_0"].x, capacitor_pos.y - flip * (150 + self.a / 2 + self.b)
            ),
        )
        _, tee_refpoints = self.insert_cell(cell_cross, cross_trans1)

        tee_port_1 = "port_right" if res_idx % 2 else "port_left"
        tee_port_2 = "port_left" if res_idx % 2 else "port_right"
        wg_1 = self._produce_waveguide(
            [
                capacitor_pos,
                pya.DPoint(capacitor_pos.x, tee_refpoints[tee_port_1].y),

            ]
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
                pya.DPoint(tee_refpoints["port_bottom"].x - 300, tee_refpoints["port_bottom"].y),
                pya.DPoint(tee_refpoints["port_bottom"].x - 500, tee_refpoints["port_bottom"].y),
                pya.DPoint(tee_refpoints["port_bottom"].x - 500, tee_refpoints["port_bottom"].y - flip * 100),
            ],
            turn_radius=turn_radius,
        )
        meander_length = total_length - wg_1.length() - cross_length - wg_2.length() - wg_3.length()
        num_meanders = _get_num_meanders(meander_length, turn_radius, w)
        self.insert_cell(
            Meander,
            start_point=pya.DPoint(tee_refpoints["port_bottom"].x - 500, tee_refpoints["port_bottom"].y - flip * 100),
            end_point=pya.DPoint(tee_refpoints["port_bottom"].x - 500, tee_refpoints["port_bottom"].y - flip * 750),
            length=meander_length,
            meanders=num_meanders,
            r=turn_radius,
        )

    def _produce_half_readout_resonator(self, cplr, cplr_dtrans, i):
        """Assemble a half-wave readout resonator; includes optional x-offset.

        A meander is inserted to complete the target length. The resulting
        resonator length is recorded in ``_readout_structure_info``.
        """
        total_length = float(self.readout_res_lengths[i])
        turn_radius = 50
        cplr_pos = self.get_refpoints(cplr, cplr_dtrans)["port_a"]

        meander_end = self._qubit_refpoints[i]["port_0"]

        if self.x_offset[i] != 0:
            flip = -1 if i % 2 else 1
            wg_1 = self._produce_waveguide(
                [
                    cplr_pos,
                    pya.DPoint(cplr_pos.x, cplr_pos.y - flip * 200),
                    pya.DPoint(cplr_pos.x - self.x_offset[i], cplr_pos.y - flip * 200),
                    pya.DPoint(cplr_pos.x - self.x_offset[i], cplr_pos.y - flip * 300),

                ]
            )
            meander_start = pya.DPoint(cplr_pos.x - self.x_offset[i], cplr_pos.y - flip * 300)
            w = 1000
            num_meanders = _get_num_meanders(total_length - wg_1.length(), turn_radius, w)
            self.insert_cell(
                Meander,
                start_point=meander_start,
                end_point=meander_end,
                length=total_length - wg_1.length(),
                meanders=num_meanders,
                r=turn_radius,
            )
        else:
            w = 1000
            num_meanders = _get_num_meanders(total_length, turn_radius, w)
            self.insert_cell(
                Meander,
                start_point=cplr_pos,
                end_point=meander_end,
                length=total_length,
                meanders=num_meanders,
                r=turn_radius,
            )
        self._readout_structure_info["readout_res_lengths"].append(total_length)

    def _produce_feedline_resonator(self):
        """Create feedline routing using a chain of tees and a meandered feedline.

        This follows the original assembly: create tees, a feedline coupling
        capacitor, then route the coplanar waveguides to connect launchers.
        """
        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a
            ),
        )
        tee_refpoints = []
        tee_rotations = [0, 2, 0, 2]
        for i in range(4):
            cross_trans = pya.DTrans(
                tee_rotations[i], False, self._qubit_refpoints[i]["port_0"].x + self.x_offset[i], 3750
            )
            inst_cross, _ = self.insert_cell(cell_cross, cross_trans)
            tee_refpoints.append(self.get_refpoints(cell_cross, inst_cross.dtrans))
        self._readout_structure_info["tees"] = [self.a, self.a, 2 * self.a]
        self._tee_refpoints = tee_refpoints

        feedline_tee_trans = pya.DTrans(0, False, 6200, 3750)
        feedline_tee, feedline_tee_refp = self.insert_cell(cell_cross, feedline_tee_trans)

        if isinstance(self.feedline_capacitor_n_fingers, tuple):
            self.feedline_capacitor_n_fingers = self.feedline_capacitor_n_fingers[0]
        cplr_params = cap_params(
            self.feedline_capacitor_n_fingers, float(self.l_fingers[0]), "smooth", finger_gap=self.b
        )
        cplr = self.add_element(**cplr_params)
        cplr_refpoints_rel = self.get_refpoints(cplr)
        cplr_dtrans = pya.DTrans(0, False, self.launchers["W"][0].x + 200, self.launchers["W"][0].y)
        _, cplr_refpoints = self.insert_cell(cplr, cplr_dtrans, rec_levels=None)

        self._produce_waveguide([
            self.launchers["W"][0],
            cplr_refpoints["port_a"],
        ], object="feedline")
        self._produce_waveguide([
            cplr_refpoints["port_b"],
            pya.DPoint(cplr_refpoints["port_b"].x + 200, cplr_refpoints["port_b"].y),
            pya.DPoint(cplr_refpoints["port_b"].x + 200, 3750 + self.feedline_meander_height),
            pya.DPoint(cplr_refpoints["port_b"].x + 400, 3750 + self.feedline_meander_height),
            pya.DPoint(cplr_refpoints["port_b"].x + 400, 3750),
            pya.DPoint(cplr_refpoints["port_b"].x + 600, 3750),
            pya.DPoint(cplr_refpoints["port_b"].x + 600, 3750 + self.feedline_meander_height),
            pya.DPoint(cplr_refpoints["port_b"].x + 800, 3750 + self.feedline_meander_height),
            pya.DPoint(cplr_refpoints["port_b"].x + 800, 3750),
            tee_refpoints[0]["port_left"],
        ], object="feedline")

        self._produce_waveguide([tee_refpoints[0]["port_right"], tee_refpoints[1]["port_right"]], object="feedline")
        self._produce_waveguide([tee_refpoints[1]["port_left"], tee_refpoints[2]["port_left"]], object="feedline")
        self._produce_waveguide([tee_refpoints[2]["port_right"], tee_refpoints[3]["port_right"]], object="feedline")
        self._produce_waveguide([tee_refpoints[3]["port_left"], feedline_tee_refp["port_left"]], object="feedline")
        self._produce_waveguide([feedline_tee_refp["port_right"], self.launchers["E"][0]], object="feedline")
        if self.tail_variant == "1":
            self._produce_waveguide([
                feedline_tee_refp["port_bottom"],
                pya.DPoint(feedline_tee_refp["port_bottom"].x, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 100, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 200, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 200, feedline_tee_refp["port_bottom"].y - 1150),
            ], object="feedline")
        elif self.tail_variant == "2":
            self._produce_waveguide([
                feedline_tee_refp["port_bottom"],
                pya.DPoint(feedline_tee_refp["port_bottom"].x, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 100, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 200, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 200, feedline_tee_refp["port_bottom"].y - 950),
            ], object="feedline")

    def _produce_readout_resonators(self):
        """Create four readout resonators connected to the tees and qubits.

        For each resonator a coupler is created and placed; depending on
        ``resonator_type`` either quarter- or half-wave building routines are
        invoked.
        """
        for i in range(4):
            cplr_params = cap_params(
                float(self.n_fingers[i]), float(self.l_fingers[i]), self.type_coupler[i], finger_gap=self.b
            )
            cplr = self.add_element(**cplr_params)
            cplr_refpoints_rel = self.get_refpoints(cplr)
            if i % 2 == 0:
                cplr_pos = (
                    self._tee_refpoints[i]["port_bottom"] - pya.DTrans.R90 * cplr_refpoints_rel["port_b"]
                )
            else:
                cplr_pos = (
                    self._tee_refpoints[i]["port_bottom"] + pya.DTrans.R90 * cplr_refpoints_rel["port_b"]
                )
            cplr_dtrans = pya.DTrans(2 * (i % 2) + 1, False, cplr_pos.x, cplr_pos.y)
            self.insert_cell(cplr, cplr_dtrans)

            if self.resonator_type == "quarter":
                self._produce_quarter_readout_resonator(cplr, cplr_dtrans, i, self._qubit_refpoints[i])
            if self.resonator_type == "half":
                self._produce_half_readout_resonator(cplr, cplr_dtrans, i)

    def _produce_chargelines(self):
        """Create charge-drive lines connecting launchers to qubit drive ports.

        The method places two tees and routes waveguides to the qubits' drive
        ports. Tapers and short waveguide segments are inserted to match
        widths where required.
        """
        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a,
            ),
        )
        cross_length = self.b * 3 + 2 * self.a
        cross_trans_0 = pya.CplxTrans(
            1,
            0,
            False,
            pya.DPoint(self.launchers["S"][0].x, self.launchers["S"][0].y + 400 + self.b + self.a / 2),
        )
        cross_trans_1 = pya.CplxTrans(
            1,
            180,
            False,
            pya.DPoint(self.launchers["N"][0].x, self.launchers["N"][0].y - 400 + self.b + self.a / 2),
        )
        _, tee_refpoints_0 = self.insert_cell(cell_cross, cross_trans_0)
        _, tee_refpoints_1 = self.insert_cell(cell_cross, cross_trans_1)
        self._produce_waveguide([
            self.launchers["S"][0],
            tee_refpoints_0["port_bottom"],
        ], a=self.a, b=self.b)
        self._produce_waveguide([
            self.launchers["N"][0],
            tee_refpoints_1["port_bottom"],
        ], a=self.a, b=self.b)

        self._produce_waveguide([
            tee_refpoints_0["port_left"],
            pya.DPoint(self._qubit_refpoints[0]["port_drive"].x, tee_refpoints_0["port_left"].y),
            pya.DPoint(self._qubit_refpoints[0]["port_drive"].x, tee_refpoints_0["port_left"].y + 100),
            pya.DPoint(self._qubit_refpoints[0]["port_drive"].x, tee_refpoints_0["port_left"].y + 400),
        ], a=self.a, b=self.b)

        self._produce_waveguide([
            tee_refpoints_1["port_right"],
            pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1["port_left"].y),
            pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1["port_left"].y - 100),
            pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1["port_left"].y - 400),
        ], a=self.a, b=self.b)

        taper_cell, _ = WaveguideCoplanarTaper.create_with_refpoints(
            self.layout,
            self.LIBRARY_NAME,
            a=self.a,
            b=self.b,
            a2=self.a / 3,
            b2=self.b / 3,
            taper_length=80,
        )
        self.insert_cell(
            taper_cell,
            pya.CplxTrans(1, 90, False, pya.DPoint(self._qubit_refpoints[0]["port_drive"].x, tee_refpoints_0["port_left"].y + 400)),
        )
        self.insert_cell(
            taper_cell,
            pya.CplxTrans(1, 90 + 180, False, pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1["port_right"].y - 400)),
        )

        self._produce_waveguide([
            pya.DPoint(self._qubit_refpoints[0]["port_drive"].x, tee_refpoints_0["port_left"].y + 480),
            self._qubit_refpoints[0]["port_drive"],
        ], a=self.a / 3, b=self.b / 3, term2=self.b)
        self._produce_waveguide([
            pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1["port_right"].y - 480),
            self._qubit_refpoints[1]["port_drive"],
        ], a=self.a / 3, b=self.b / 3, term2=self.b)

        self._produce_waveguide([
            tee_refpoints_0["port_right"],
            pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0["port_left"].y),
            pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0["port_left"].y + 100),
            pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0["port_left"].y + 400),
        ], a=self.a, b=self.b)

        self._produce_waveguide([
            tee_refpoints_1["port_left"],
            pya.DPoint(self._qubit_refpoints[3]["port_drive"].x, tee_refpoints_1["port_left"].y),
            pya.DPoint(self._qubit_refpoints[3]["port_drive"].x, tee_refpoints_1["port_left"].y - 100),
            pya.DPoint(self._qubit_refpoints[3]["port_drive"].x, tee_refpoints_1["port_left"].y - 400),
        ], a=self.a, b=self.b)
        self.insert_cell(
            taper_cell,
            pya.CplxTrans(1, 90, False, pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0["port_left"].y + 400)),
        )
        self.insert_cell(
            taper_cell,
            pya.CplxTrans(1, 90 + 180, False, pya.DPoint(self._qubit_refpoints[3]["port_drive"].x, tee_refpoints_1["port_right"].y - 400)),
        )
        self._produce_waveguide([
            pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0["port_left"].y + 480),
            self._qubit_refpoints[2]["port_drive"],
        ], a=self.a / 3, b=self.b / 3, term2=self.b)
        self._produce_waveguide([
            pya.DPoint(self._qubit_refpoints[3]["port_drive"].x, tee_refpoints_1["port_right"].y - 480),
            self._qubit_refpoints[3]["port_drive"],
        ], a=self.a / 3, b=self.b / 3, term2=self.b)

    def _produce_chargelines_v2(self):
        """Alternate charge-line routing variant used when requested by params."""
        tee_refpoints_0_placeholder = pya.DPoint(self.launchers["S"][0].x, self.launchers["S"][0].y + 400 + self.b + self.a / 2)
        tee_refpoints_1_placeholder = pya.DPoint(self.launchers["N"][0].x, self.launchers["N"][0].y - 400 + self.b + self.a / 2)
        self._produce_waveguide([
            self.launchers["S"][0],
            tee_refpoints_0_placeholder,
            pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0_placeholder.y),
            pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0_placeholder.y + 100),
            pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0_placeholder.y + 400),
        ], a=self.a, b=self.b)

        self._produce_waveguide([
            self.launchers["N"][0],
            tee_refpoints_1_placeholder,
            pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1_placeholder.y),
            pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1_placeholder.y - 100),
            pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1_placeholder.y - 400),
        ], a=self.a, b=self.b)

        taper_cell, _ = WaveguideCoplanarTaper.create_with_refpoints(
            self.layout,
            self.LIBRARY_NAME,
            a=self.a,
            b=self.b,
            a2=self.a / 3,
            b2=self.b / 3,
            taper_length=80,
        )

        self.insert_cell(
            taper_cell,
            pya.CplxTrans(1, 90 + 180, False, pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1_placeholder.y - 400)),
        )

        self._produce_waveguide([
            pya.DPoint(self._qubit_refpoints[1]["port_drive"].x, tee_refpoints_1_placeholder.y - 480),
            self._qubit_refpoints[1]["port_drive"],
        ], a=self.a / 3, b=self.b / 3, term2=self.b)

        self.insert_cell(
            taper_cell,
            pya.CplxTrans(1, 90, False, pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0_placeholder.y + 400)),
        )

        self._produce_waveguide([
            pya.DPoint(self._qubit_refpoints[2]["port_drive"].x, tee_refpoints_0_placeholder.y + 480),
            self._qubit_refpoints[2]["port_drive"],
        ], a=self.a / 3, b=self.b / 3, term2=self.b)

    def _produce_feedline(self, x_distance):
        """Create a simplified feedline (non-meander resonator) variant."""
        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a
            ),
        )
        tee_refpoints = []
        tee_rotations = [0, 2, 0, 2]
        for i in range(4):
            cross_trans = pya.DTrans(
                tee_rotations[i], False, self._qubit_refpoints[i]["port_0"].x + float(self.x_offset[i]), 3.75e3
            )
            inst_cross, _ = self.insert_cell(cell_cross, cross_trans)
            tee_refpoints.append(self.get_refpoints(cell_cross, inst_cross.dtrans))
        self._readout_structure_info["tees"] = [self.a, self.a, 2 * self.a]
        self._tee_refpoints = tee_refpoints
        wg0 = self._produce_waveguide([
            self.launchers["W"][0],
            pya.DPoint(self.launchers["W"][0].x + x_distance, self.launchers["W"][0].y),
            tee_refpoints[0]["port_left"],
        ])
        self._readout_structure_info["feedline"].append(wg0.length())

        wg1 = self._produce_waveguide([tee_refpoints[0]["port_right"], tee_refpoints[1]["port_right"]])
        self._readout_structure_info["feedline"].append(wg1.length())

        wg2 = self._produce_waveguide([tee_refpoints[1]["port_left"], tee_refpoints[2]["port_left"]])
        self._readout_structure_info["feedline"].append(wg2.length())

        wg3 = self._produce_waveguide([tee_refpoints[2]["port_right"], tee_refpoints[3]["port_right"]])
        self._readout_structure_info["feedline"].append(wg3.length())

        wg4 = self._produce_waveguide([
            tee_refpoints[3]["port_left"],
            pya.DPoint(self.launchers["E"][0].x - x_distance, self.launchers["E"][0].y),
            self.launchers["E"][0],
        ])
        self._readout_structure_info["feedline"].append(wg4.length())
