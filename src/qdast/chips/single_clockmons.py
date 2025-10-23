"""
PCell implementation for SingleClockmons chip (4 Clockmon qubits around a feedline).
"""
from typing import Tuple
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


def _get_num_meanders(
    meander_length: float, turn_radius: float, meander_min_width: float
) -> int:
    """Return the integer number of meanders required for a given length.

    Uses the same formula as the original implementation: available straight
    length divided by the effective per-meander width (including arc length
    from turns).
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
class SingleClockmons(QDASTChip):
    """PCell declaration for a SingleClockmons chip.

    The device places four Clockmon qubits coupled to a common feedline.
    Half of the qubits are above and half are below the feedline. Default
    parameters and placement mirror the original source.
    """

    readout_res_lengths = Param(
        pdt.TypeList,
        "Readout resonator lengths (four resonators)",
        [7000, 7100, 7200, 7300],
    )
    n_fingers = Param(
        pdt.TypeList,
        "Number of fingers for readout resonator couplers",
        [1.2, 1.3, 1.4, 1.5],
    )
    feedline_capacitor_n_fingers = Param(
        pdt.TypeDouble, "Number of fingers for feedline resonator coupler", 3.1
    )

    type_coupler = Param(
        pdt.TypeList,
        "Coupler type for test resonator couplers",
        ["smooth", "smooth", "smooth", "smooth"],
    )
    coupler_widths = Param(pdt.TypeList, "Qubit coupler width", [150, 150, 150, 150])

    _readout_structure_info = {"feedline": [], "tees": [], "readout_res_lengths": []}

    l_fingers = Param(
        pdt.TypeList,
        "Length of fingers for readout resonator couplers",
        [30, 30, 30, 30],
    )
    with_feedline_resonator = Param(pdt.TypeBoolean, "Enable feedline resonator", True)

    def build(self) -> None:
        """Top-level assembly called by the PCell framework.

        Places launchers, qubits and then assembles the chosen feedline
        variant and the readout resonators.
        """
        self.launchers = self.produce_launchers(
            "SquareNSWE_5x5", launcher_assignments={4: "FL-IN", 2: "FL-OUT"}
        )
        self.qubits_refpoints = self._produce_qubits()

        feedline_x_distance = 200
        if self.with_feedline_resonator:
            self._produce_feedline_resonator(feedline_x_distance)
        else:
            self._produce_feedline(feedline_x_distance)

        self._produce_readout_resonators()

    def _produce_waveguide(self, path, term2: int = 0, turn_radius=None):
        """Create a WaveguideCoplanar element following ``path`` and insert it.

        The created waveguide length is appended to ``_readout_structure_info['feedline']``.
        """
        if turn_radius is None:
            turn_radius = self.r
        waveguide = self.add_element(
            WaveguideCoplanar, path=pya.DPath(path, 1), r=turn_radius, term2=term2
        )
        self.insert_cell(waveguide)
        self._readout_structure_info["feedline"].append(waveguide.length())
        return waveguide

    def _produce_qubit(
        self,
        coupler_width: float,
        center_x: float,
        center_y: float,
        rotation: int,
        name: str = None,
    ):
        """Instantiate a Clockmon qubit and insert it at the requested location.

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
            clock_diameter=95,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered",
            sim_tool=self.sim_tool,
            with_squid=self.with_squid,
            pad_width=6,
            taper_width=95 / 7,
        )
        qubit_trans = pya.DTrans(rotation, False, center_x, center_y)
        _, refpoints_abs = self.insert_cell(qubit, qubit_trans, name, rec_levels=None)
        return refpoints_abs

    def _produce_qubits(self) -> Tuple:
        """Create four qubits and store their refpoints.

        Returns a tuple of the four qubit refpoints (qb0..qb3).
        """
        qubit_spacing_x = 850
        qubit_spacing_x_alt = 500
        qubit_spacing_y = 1200
        qubits_center_x = 2.5e3 if self.with_feedline_resonator else 2.2e3

        y_a = 3.5e3 + qubit_spacing_y / 2
        y_b = 1.5e3 - qubit_spacing_y / 2

        qb0_refpoints = self._produce_qubit(
            float(self.coupler_widths[0]),
            qubits_center_x - qubit_spacing_x,
            y_b,
            0,
            "qb_0",
        )
        qb1_refpoints = self._produce_qubit(
            float(self.coupler_widths[1]),
            qubits_center_x - qubit_spacing_x + qubit_spacing_x_alt + 100,
            y_a,
            2,
            "qb_1",
        )
        qb2_refpoints = self._produce_qubit(
            float(self.coupler_widths[2]),
            qubits_center_x + qubit_spacing_x - 50,
            y_b,
            0,
            "qb_2",
        )
        qb3_refpoints = self._produce_qubit(
            float(self.coupler_widths[3]),
            qubits_center_x + qubit_spacing_x + qubit_spacing_x_alt - 50,
            y_a,
            2,
            "qb_3",
        )

        self._qubit_x_coords = [
            qubits_center_x - qubit_spacing_x,
            qubits_center_x - qubit_spacing_x + qubit_spacing_x_alt,
            qubits_center_x + qubit_spacing_x,
            qubits_center_x + qubit_spacing_x + qubit_spacing_x_alt,
        ]
        self._qubit_refpoints = [
            qb0_refpoints,
            qb1_refpoints,
            qb2_refpoints,
            qb3_refpoints,
        ]
        return qb0_refpoints, qb1_refpoints, qb2_refpoints, qb3_refpoints

    def _produce_readout_resonator(self, capacitor, capacitor_dtrans, res_idx: int):
        """Assemble a readout resonator by inserting a Meander between the capacitor and qubit port.

        The total length is taken from ``readout_res_lengths[res_idx]``. The method
        records the chosen length in ``_readout_structure_info``.
        """
        total_length = float(self.readout_res_lengths[res_idx])
        turn_radius = 50

        meander_start = self.get_refpoints(capacitor, capacitor_dtrans)["port_a"]
        meander_end = self._qubit_refpoints[res_idx]["port_0"]

        w = 1000
        num_meanders = _get_num_meanders(total_length, turn_radius, w)
        self.insert_cell(
            Meander,
            start_point=meander_start,
            end_point=meander_end,
            length=total_length,
            meanders=num_meanders,
            r=turn_radius,
        )
        self._readout_structure_info["readout_res_lengths"].append(total_length)

    def _produce_feedline(self, x_distance: float):
        """Produce a simple feedline connecting FL-IN to FL-OUT via four tees."""
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
                tee_rotations[i], False, self._qubit_refpoints[i]["port_0"].x, 2.5e3
            )
            inst_cross, _ = self.insert_cell(cell_cross, cross_trans)
            tee_refpoints.append(self.get_refpoints(cell_cross, inst_cross.dtrans))
        self._readout_structure_info["tees"] = [self.a, self.a, 2 * self.a]
        self._tee_refpoints = tee_refpoints

        self._produce_waveguide(
            [
                self.launchers["FL-IN"][0],
                pya.DPoint(
                    self.launchers["FL-IN"][0].x + x_distance,
                    self.launchers["FL-IN"][0].y,
                ),
                tee_refpoints[0]["port_left"],
            ]
        )
        self._produce_waveguide(
            [tee_refpoints[0]["port_right"], tee_refpoints[1]["port_right"]]
        )
        self._produce_waveguide(
            [tee_refpoints[1]["port_left"], tee_refpoints[2]["port_left"]]
        )
        self._produce_waveguide(
            [tee_refpoints[2]["port_right"], tee_refpoints[3]["port_right"]]
        )
        self._produce_waveguide(
            [
                tee_refpoints[3]["port_left"],
                pya.DPoint(
                    self.launchers["FL-OUT"][0].x - x_distance,
                    self.launchers["FL-OUT"][0].y,
                ),
                self.launchers["FL-OUT"][0],
            ]
        )

    def _produce_feedline_resonator(self, x_distance: float):
        """Produce a feedline variant which includes a coupling capacitor and special routing.

        The function keeps the same geometry and routing choices as the original
        implementation; it also records feedline segment lengths in
        ``_readout_structure_info``.
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
                tee_rotations[i], False, self._qubit_refpoints[i]["port_0"].x, 2.5e3
            )
            inst_cross, _ = self.insert_cell(cell_cross, cross_trans)
            tee_refpoints.append(self.get_refpoints(cell_cross, inst_cross.dtrans))
        self._readout_structure_info["tees"] = [self.a, self.a, 2 * self.a]
        self._tee_refpoints = tee_refpoints

        feedline_tee_trans = pya.DTrans(0, False, 4100, 2.5e3)
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
            0, False, self.launchers["FL-IN"][0].x + 200, self.launchers["FL-IN"][0].y
        )
        _, cplr_refpoints = self.insert_cell(cplr, cplr_dtrans, rec_levels=None)

        self._produce_waveguide([self.launchers["FL-IN"][0], cplr_refpoints["port_a"]])
        self._produce_waveguide(
            [
                cplr_refpoints["port_b"],
                pya.DPoint(cplr_refpoints["port_b"].x + 80, cplr_refpoints["port_b"].y),
                pya.DPoint(cplr_refpoints["port_b"].x + 80, 4530),
                pya.DPoint(cplr_refpoints["port_b"].x + 180, 4530),
                pya.DPoint(cplr_refpoints["port_b"].x + 180, 2500),
                pya.DPoint(cplr_refpoints["port_b"].x + 300, 2500),
                pya.DPoint(cplr_refpoints["port_b"].x + 300, 4530),
                pya.DPoint(cplr_refpoints["port_b"].x + 400, 4530),
                pya.DPoint(cplr_refpoints["port_b"].x + 400, 2500),
                tee_refpoints[0]["port_left"],
            ],
            turn_radius=30,
        )

        self._produce_waveguide(
            [tee_refpoints[0]["port_right"], tee_refpoints[1]["port_right"]]
        )
        self._produce_waveguide(
            [tee_refpoints[1]["port_left"], tee_refpoints[2]["port_left"]]
        )
        self._produce_waveguide(
            [tee_refpoints[2]["port_right"], tee_refpoints[3]["port_right"]]
        )
        self._produce_waveguide(
            [tee_refpoints[3]["port_left"], feedline_tee_refp["port_left"]]
        )
        self._produce_waveguide(
            [feedline_tee_refp["port_right"], self.launchers["FL-OUT"][0]]
        )
        self._produce_waveguide(
            [
                feedline_tee_refp["port_bottom"],
                pya.DPoint(
                    feedline_tee_refp["port_bottom"].x,
                    feedline_tee_refp["port_bottom"].y - 500,
                ),
                pya.DPoint(
                    feedline_tee_refp["port_bottom"].x + 100,
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
            ]
        )

    def _produce_readout_resonators(self) -> None:
        """Create readout resonators for each coupler/tee pair and insert them."""
        for i in range(4):
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

            self._produce_readout_resonator(cplr, cplr_dtrans, i)

    # Utility to export readout structure info
    # def get_readout_structure_info(self):
    #     if self.with_feedline_resonator:
    #         self._readout_structure_info["readout_structure"] = ["feedline_resonator"]
    #     else:
    #         self._readout_structure_info["readout_structure"] = ["simple"]
    #     df = pd.DataFrame.from_dict(self._readout_structure_info, orient = "index")
    #     df.to_csv("single_clockmons_readout_structure.csv", mode='w+')
