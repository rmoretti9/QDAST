"""
PCell for a two-Clockmon chip used to study capacitive coupling and CR-like
interactions.

This module is a light refactor of the original source: docstrings were
added and helper functions slightly reorganized while preserving the
original logic, parameter names and layout semantics.
"""

from math import pi
import os
from typing import Tuple, Dict

import pandas as pd
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.elements.meander import Meander
from qdast.qubits.clockmon import Clockmon
from kqcircuits.elements.waveguide_coplanar import WaveguideCoplanar
from kqcircuits.elements.waveguide_coplanar_taper import WaveguideCoplanarTaper

from kqcircuits.elements.waveguide_coplanar_splitter import (
    WaveguideCoplanarSplitter,
    t_cross_parameters,
)
from kqcircuits.pya_resolver import pya
from kqcircuits.util.coupler_lib import cap_params
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
from kqcircuits.junctions.junction import Junction


def _get_num_meanders(meander_length, turn_radius, meander_min_width):
    """Return the integer number of meanders required for a given length.

    The calculation mirrors the original implementation: divide the
    available straight length by the effective width per meander including
    turning arcs.
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
class TwoClockmons(QDASTChip):
    """PCell declaration for a two-Clockmons chip.

    The device places two Clockmon qubits capacitively coupled to a bus and
    includes readout resonators and drive lines. Parameters and behavior are
    unchanged from the original implementation; docstrings added for
    readability.
    """

    readout_res_lengths = Param(pdt.TypeList, "Readout resonator lengths", [4000, 4100])
    n_fingers = Param(
        pdt.TypeList, "Number of fingers for readout resonator couplers", [1.2, 1.3]
    )
    feedline_capacitor_n_fingers = Param(
        pdt.TypeDouble, "Number of fingers for feedline resonator coupler", 3.1
    )

    type_coupler = Param(
        pdt.TypeList, "Coupler type for test resonator couplers", ["smooth", "smooth"]
    )
    readout_coupler_widths = Param(pdt.TypeList, "Readout coupler width", [150, 150])
    qubit_coupler_widths = Param(pdt.TypeList, "Qubit coupler width", [70, 70])

    l_fingers = Param(
        pdt.TypeList, "Length of fingers for readout resonator couplers", [30, 30]
    )
    with_feedline_resonator = Param(pdt.TypeBoolean, "Enable feedline resonator", False)
    coupler_length = Param(pdt.TypeDouble, "Bus coupler length", 8000)

    def build(self):
        """Top-level assembly called by the PCell framework.

        Produces launchers, qubits, feedline, readout resonators, bus
        coupler and charge-drive lines.
        """
        self._readout_structure_info: Dict = {}
        self.launchers = self.produce_launchers(
            "SquareNSWE_5x5",
            launcher_assignments={4: "FL-IN", 1: "FL-OUT", 3: "DL-0", 2: "DL-1"},
        )
        self.qubits_refpoints = self._produce_qubits()

        self._produce_feedline()
        self._produce_readout_resonators()
        self._produce_coupler()
        self._produce_chargelines()

    def _produce_waveguide(self, path, term2=0, turn_radius=None, a=None, b=None):
        """Add a WaveguideCoplanar element that follows ``path`` and insert it.

        Returns the created waveguide element (so callers can query .length()).
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
        return waveguide

    def _produce_qubit(
        self,
        qubit_id: int,
        center_x: float,
        center_y: float,
        rotation: int,
        name: str = None,
    ):
        """Instantiate a Clockmon qubit and insert it into the layout.

        Returns the absolute reference points created by the insertion.
        """
        rr_width = float(self.readout_coupler_widths[qubit_id])
        qb_width = float(self.qubit_coupler_widths[qubit_id])
        coupler_widths = (
            [rr_width, 0, qb_width, 0, 0, 0]
            if qubit_id == 0
            else [rr_width, 0, 0, qb_width, 0, 0]
        )
        coupler_offsets = (
            [255, 0, 288, 0, 0, 0] if qubit_id == 0 else [255, 0, 0, 250, 0, 0]
        )

        qubit = self.add_element(
            Clockmon,
            ground_gap=[630, 610],
            a=10,
            b=6,
            island_extent=[535, 200],
            coupler_widths=coupler_widths,
            coupler_offsets=coupler_offsets,
            island_to_island_distance=50,
            clock_diameter=95,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered",
            sim_tool=self.sim_tool,
            with_squid=self.with_squid,
            pad_width=6,
            taper_width=95 / 7,
            drive_position=[-400, -550] if qubit_id == 0 else [450, -200],
        )
        qubit_trans = pya.CplxTrans(1, rotation, False, center_x, center_y)
        _, refpoints_abs = self.insert_cell(qubit, qubit_trans, name, rec_levels=None)
        return refpoints_abs

    def _produce_qubits(self) -> Tuple[Dict, Dict]:
        """Create two Clockmon qubits at predefined positions and return their refpoints."""
        qb0_refpoints = self._produce_qubit(0, 1600, 1600, 0, "qb_0")
        qb1_refpoints = self._produce_qubit(1, 3200, 3000, 0, "qb_1")
        self._qubit_refpoints = [qb0_refpoints, qb1_refpoints]
        return qb0_refpoints, qb1_refpoints

    def _produce_readout_resonator(self, capacitor, capacitor_dtrans, res_idx: int):
        """Assemble a readout resonator for the specified index and attach it to a tee.

        A meander section is added to reach the requested total length and
        intermediate segment lengths are recorded in ``_readout_structure_info``.
        """
        total_length = float(self.readout_res_lengths[res_idx])
        turn_radius = 50
        w = 800

        wg_start = self.get_refpoints(capacitor, capacitor_dtrans)["port_a"]
        qubit_port = self._qubit_refpoints[res_idx]["port_0"]

        if res_idx == 0:
            cell_cross = self.add_element(
                WaveguideCoplanarSplitter,
                **t_cross_parameters(
                    a=self.a,
                    b=self.b,
                    a2=self.a,
                    b2=self.b,
                    length_extra_side=2 * self.a,
                ),
            )
            cross_length = self.b * 3 + 2 * self.a
            wg_1 = self._produce_waveguide(
                [
                    wg_start,
                    pya.DPoint(
                        self.qubits_refpoints[res_idx]["port_0"].x - self.a - 300,
                        wg_start.y,
                    ),
                ]
            )
            cross_trans1 = pya.CplxTrans(
                1,
                0,
                False,
                pya.DPoint(
                    self.qubits_refpoints[res_idx]["port_0"].x - 300, wg_start.y
                ),
            )
            _, tee_refpoints = self.insert_cell(cell_cross, cross_trans1)
            wg_2 = self._produce_waveguide(
                [
                    tee_refpoints["port_bottom"],
                    pya.DPoint(qubit_port.x - 300, qubit_port.y + 300),
                    pya.DPoint(qubit_port.x, qubit_port.y + 300),
                    pya.DPoint(qubit_port.x, qubit_port.y),
                ]
            )
            wg_3 = self._produce_waveguide(
                [
                    tee_refpoints["port_right"],
                    pya.DPoint(
                        tee_refpoints["port_right"].x + 300,
                        tee_refpoints["port_right"].y,
                    ),
                    pya.DPoint(
                        tee_refpoints["port_right"].x + 350,
                        tee_refpoints["port_right"].y + 50,
                    ),
                ]
            )
            meander_length = (
                total_length
                - wg_1.length()
                - cross_length
                - wg_2.length()
                - wg_3.length()
            )
            num_meanders = _get_num_meanders(meander_length, turn_radius, w)
            self.insert_cell(
                Meander,
                start_point=pya.DPoint(
                    tee_refpoints["port_right"].x + 350,
                    tee_refpoints["port_right"].y + 50,
                ),
                end_point=pya.DPoint(
                    tee_refpoints["port_right"].x + 700,
                    tee_refpoints["port_right"].y + 400,
                ),
                length=meander_length,
                meanders=num_meanders,
                r=turn_radius,
            )

        if res_idx == 1:
            cell_cross = self.add_element(
                WaveguideCoplanarSplitter,
                **t_cross_parameters(
                    a=self.a,
                    b=self.b,
                    a2=self.a,
                    b2=self.b,
                    length_extra_side=2 * self.a,
                    angles=[90, 270, 0],
                ),
            )
            cross_length = self.b * 3 + 2 * self.a
            wg_1 = self._produce_waveguide(
                [wg_start, pya.DPoint(wg_start.x, wg_start.y - 200)]
            )
            cross_trans1 = pya.CplxTrans(
                1, 0, False, pya.DPoint(wg_start.x, wg_start.y - 200 - self.a / 2)
            )
            _, tee_refpoints = self.insert_cell(cell_cross, cross_trans1)
            wg_2 = self._produce_waveguide(
                [
                    tee_refpoints["port_left"],
                    pya.DPoint(
                        tee_refpoints["port_left"].x, tee_refpoints["port_left"].y - 200
                    ),
                    pya.DPoint(qubit_port.x, tee_refpoints["port_left"].y - 200),
                    qubit_port,
                ]
            )
            wg_3 = self._produce_waveguide(
                [
                    tee_refpoints["port_bottom"],
                    pya.DPoint(
                        tee_refpoints["port_bottom"].x + 800,
                        tee_refpoints["port_bottom"].y,
                    ),
                ]
            )
            meander_length = (
                total_length
                - wg_1.length()
                - cross_length
                - wg_2.length()
                - wg_3.length()
            )
            num_meanders = _get_num_meanders(meander_length, turn_radius, w)
            self.insert_cell(
                Meander,
                start_point=pya.DPoint(
                    tee_refpoints["port_bottom"].x + 800, tee_refpoints["port_bottom"].y
                ),
                end_point=pya.DPoint(
                    tee_refpoints["port_bottom"].x + 850 + 550,
                    tee_refpoints["port_bottom"].y,
                ),
                length=meander_length,
                meanders=num_meanders,
                r=turn_radius,
            )

        self._readout_structure_info["readout_res_" + str(res_idx + 1)] = (
            wg_1.length(),
            wg_2.length(),
            wg_3.length(),
            meander_length,
        )

    def _produce_feedline(self):
        """Create a simple feedline routing between two tee locations."""
        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a
            ),
        )

        cross_trans1 = pya.CplxTrans(1, 90, False, 900, 2.7e3)
        _, tee_refpoints1 = self.insert_cell(cell_cross, cross_trans1)

        cross_trans2 = pya.CplxTrans(
            1, 0, False, 2.0e3, self.launchers["FL-OUT"][0].y - 200
        )
        _, tee_refpoints2 = self.insert_cell(cell_cross, cross_trans2)

        tee_refpoints = [tee_refpoints1, tee_refpoints2]

        self._readout_structure_info["tees"] = [self.a, self.a, 2 * self.a]
        self._tee_refpoints = tee_refpoints

        fl_1 = self._produce_waveguide(
            [
                self.launchers["FL-IN"][0],
                pya.DPoint(tee_refpoints[0]["port_left"].x, 2.5e3),
                tee_refpoints[0]["port_left"],
            ]
        )
        fl_2 = self._produce_waveguide(
            [
                tee_refpoints[0]["port_right"],
                pya.DPoint(
                    tee_refpoints[0]["port_right"].x, tee_refpoints[1]["port_left"].y
                ),
                tee_refpoints[1]["port_left"],
            ]
        )
        fl_3 = self._produce_waveguide(
            [
                tee_refpoints[1]["port_right"],
                pya.DPoint(
                    self.launchers["FL-OUT"][0].x, self.launchers["FL-OUT"][0].y - 200
                ),
                self.launchers["FL-OUT"][0],
            ]
        )

        self._readout_structure_info["feedline"] = (
            fl_1.length(),
            fl_2.length(),
            fl_3.length(),
        )

    def _produce_readout_resonators(self):
        """Create two readout resonators connected to feedline tees.

        For each coupler a capacitor element is created and placed, then the
        corresponding resonator is assembled.
        """
        for i in range(2):
            cplr_params = cap_params(
                float(self.n_fingers[i]),
                float(self.l_fingers[i]),
                self.type_coupler[i],
                finger_gap=self.b,
            )
            cplr = self.add_element(**cplr_params)
            cplr_refpoints_rel = self.get_refpoints(cplr)
            smooth_capacitor_height = abs(
                cplr_refpoints_rel["base"].x - cplr_refpoints_rel["port_b"].x
            )
            if i == 0:
                cplr_pos = self._tee_refpoints[i]["port_bottom"]
                cplr_dtrans = pya.DCplxTrans(
                    1, 180, False, cplr_pos.x + smooth_capacitor_height, cplr_pos.y
                )
            elif i == 1:
                cplr_pos = self._tee_refpoints[i]["port_bottom"]
                cplr_dtrans = pya.DCplxTrans(
                    1, 90, False, cplr_pos.x, cplr_pos.y - smooth_capacitor_height
                )
            self.insert_cell(cplr, cplr_dtrans)
            self._produce_readout_resonator(cplr, cplr_dtrans, i)

    def _produce_coupler(self):
        """Assemble the capacitive bus coupling the two qubits.

        The method routes waveguides from each qubit to a common region and
        inserts a meander of the appropriate remaining length.
        """
        wg_qb0 = self._produce_waveguide(
            [
                self.qubits_refpoints[0]["port_2"],
                pya.DPoint(
                    self.qubits_refpoints[0]["port_2"].x + 700,
                    self.qubits_refpoints[0]["port_2"].y,
                ),
                pya.DPoint(
                    self.qubits_refpoints[0]["port_2"].x + 700,
                    self.qubits_refpoints[0]["port_2"].y - 200,
                ),
                pya.DPoint(
                    self.qubits_refpoints[0]["port_2"].x + 900,
                    self.qubits_refpoints[0]["port_2"].y - 200,
                ),
                pya.DPoint(
                    self.qubits_refpoints[1]["port_3"].x,
                    self.qubits_refpoints[0]["port_2"].y - 200,
                ),
                pya.DPoint(
                    self.qubits_refpoints[1]["port_3"].x,
                    self.qubits_refpoints[0]["port_2"].y,
                ),
            ]
        )
        wg_qb1 = self._produce_waveguide(
            [
                self.qubits_refpoints[1]["port_3"],
                pya.DPoint(
                    self.qubits_refpoints[1]["port_3"].x,
                    self.qubits_refpoints[1]["port_3"].y - 100,
                ),
            ]
        )

        meander_length = self.coupler_length - wg_qb0.length() - wg_qb1.length()
        turn_radius = 50
        w = 1000
        num_meanders = _get_num_meanders(meander_length, turn_radius, w)
        self.insert_cell(
            Meander,
            start_point=pya.DPoint(
                self.qubits_refpoints[1]["port_3"].x,
                self.qubits_refpoints[0]["port_2"].y,
            ),
            end_point=pya.DPoint(
                self.qubits_refpoints[1]["port_3"].x,
                self.qubits_refpoints[1]["port_3"].y - 100,
            ),
            length=meander_length,
            meanders=num_meanders,
            r=turn_radius,
        )

    def _produce_chargelines(self):
        """Create drive lines from launcher drive ports to each qubit drive port.

        Includes tapers to handle width transitions and returns the created
        waveguide elements for the two drive lines.
        """
        dl_0 = self._produce_waveguide(
            [
                self.launchers["DL-0"][0],
                pya.DPoint(
                    self.launchers["DL-0"][0].x, self.launchers["DL-0"][0].y + 100
                ),
                pya.DPoint(
                    self.qubits_refpoints[0]["port_drive"].x,
                    self.launchers["DL-0"][0].y + 100,
                ),
                pya.DPoint(
                    self.qubits_refpoints[0]["port_drive"].x,
                    self.launchers["DL-0"][0].y + 200,
                ),
            ]
        )

        taper_cell, taper_ref_rel = WaveguideCoplanarTaper.create_with_refpoints(
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
                    self.qubits_refpoints[0]["port_drive"].x,
                    self.launchers["DL-0"][0].y + 200,
                ),
            ),
        )
        dl_0_tapered = self._produce_waveguide(
            [taper_ref0["port_b"], self.qubits_refpoints[0]["port_drive"]],
            a=self.a / 3,
            b=self.b / 3,
            term2=self.b,
        )

        dl_1 = self._produce_waveguide(
            [
                self.launchers["DL-1"][0],
                pya.DPoint(
                    self.launchers["DL-1"][0].x - 150, self.launchers["DL-1"][0].y
                ),
                pya.DPoint(
                    self.launchers["DL-1"][0].x - 150,
                    self.qubits_refpoints[1]["port_drive"].y,
                ),
                pya.DPoint(
                    self.launchers["DL-1"][0].x - 250,
                    self.qubits_refpoints[1]["port_drive"].y,
                ),
            ]
        )
        _, taper_ref1 = self.insert_cell(
            taper_cell,
            pya.CplxTrans(
                1,
                180,
                False,
                pya.DPoint(
                    self.launchers["DL-1"][0].x - 250,
                    self.qubits_refpoints[1]["port_drive"].y,
                ),
            ),
        )
        dl_1_tapered = self._produce_waveguide(
            [taper_ref1["port_b"], self.qubits_refpoints[1]["port_drive"]],
            a=self.a / 3,
            b=self.b / 3,
            term2=self.b,
        )

        return dl_0_tapered, dl_1_tapered
