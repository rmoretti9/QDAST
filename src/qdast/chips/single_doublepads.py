from math import pi
import pandas as pd
from qdast.chips.qdast_chip import QDASTChip
from kqcircuits.elements.meander import Meander
from qdast.qubits.clockmon import Clockmon
from kqcircuits.qubits.double_pads import DoublePads
from kqcircuits.elements.waveguide_coplanar import WaveguideCoplanar
import os

from kqcircuits.elements.waveguide_coplanar_splitter import (
    WaveguideCoplanarSplitter,
    t_cross_parameters,
)
from kqcircuits.pya_resolver import pya
from kqcircuits.util.coupler_lib import cap_params
from kqcircuits.util.parameters import Param, pdt, add_parameters_from
from kqcircuits.junctions.junction import Junction
from kqcircuits.elements.launcher import Launcher

def _get_num_meanders(meander_length, turn_radius, meander_min_width):
    """Get the required number of meanders to create a meander element with the given parameters."""

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
class SingleDoublepads(QDASTChip):
    """The PCell declaration for a SingleDoublepads chip.
    """
    # name_mask = Param(pdt.TypeString, "Name of the mask", "M000")  # string '_3' will leave empty space for M000
    # name_chip = Param(pdt.TypeString, "Name of the chip", "CTest")
    # name_copy = Param(pdt.TypeString, "Name of the copy", None)
    readout_res_lengths = Param(
        pdt.TypeList,
        "Readout resonator lengths (four resonators)",
        [6900, 7000, 7100, 7200, 7300, 7400],
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
    # coupler_r = Param(pdt.TypeDouble, "Coupler rounding radius", 10, unit="μm")

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

        self.insert_cell(launcher_cell, pya.DCplxTrans(1, 90, False, 453 + 200, 10000 - 620 - 200))
        self.insert_cell(launcher_cell, pya.DCplxTrans(1, 90, False, 10000 - 453 - 200, 10000 - 615 - 200))
        _, self.refpoints_fl_out = self.insert_cell(launcher_cell, pya.DCplxTrans(1, 0, False, 10000 - 620 - 200, 5000)) #
        self.insert_cell(launcher_cell, pya.DCplxTrans(1, -90, False, 10000 - 453 - 200, + 615 + 200))
        self.insert_cell(launcher_cell, pya.DCplxTrans(1, -90, False, 453 + 200 + 7, + 615 + 200))
        _, self.refpoints_fl_in = self.insert_cell(launcher_cell, pya.DCplxTrans(1, -180, False, 620 + 200, 5000)) #

        self._produce_qubits()

        self._produce_feedline_resonator()
        self._produce_readout_resonators()
        # # self.get_readout_structure_info()

    def _produce_waveguide(self, path, term2=0, turn_radius=None):
        if turn_radius is None:
            turn_radius = self.r
        waveguide = self.add_element(
            WaveguideCoplanar,
            path=pya.DPath(path, 1),
            r=turn_radius,
            term2=term2,
        )
        self.insert_cell(waveguide)
        self._readout_structure_info["feedline"].append(waveguide.length())
        return waveguide

    def _produce_qubit(self, coupler_width, center_x, center_y, rotation, name=None):
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
            sim_tool = self.sim_tool,
            with_squid = self.with_squid,
            pad_width = 6,
            taper_width = 95/7,
            bent_section_length  = 8,
            lead_height_untapered = 4,
            lead_height_tapered = 8
        )
        qubit_trans = pya.DTrans(rotation, False, center_x, center_y)
        _, refpoints_abs = self.insert_cell(
            qubit, qubit_trans, name, rec_levels=None
        )
        return refpoints_abs

    def _produce_qubits(self):
        qubit_spacing_x = 1500  # x-distance between qubits on the same feedline side
        qubit_spacing_x_alt = 800 # x-distance between qubits on different feedline sides
        qubit_spacing_y = 5e3  # shortest y-distance between qubit centers on different sides of the feedline
        qubits_center_x = 5e3  # the x-coordinate around which qubits are centered
        # qubits above the feedline, from left to right
        y_a = 5e3 + qubit_spacing_y / 2
        y_b = 5e3 - qubit_spacing_y / 2
        qb0_refpoints = self._produce_qubit(
            float(self.coupler_widths[0]), qubits_center_x - qubit_spacing_x, y_b, 0, "qb_0"
        )
        qb1_refpoints = self._produce_qubit(
            float(self.coupler_widths[1]), qubits_center_x - qubit_spacing_x + qubit_spacing_x_alt, y_a, 2, "qb_1"
        )
        qb2_refpoints = self._produce_qubit(
            float(self.coupler_widths[2]), qubits_center_x, y_b, 0, "qb_2"
        )
        qb3_refpoints = self._produce_qubit(
            float(self.coupler_widths[3]), qubits_center_x + qubit_spacing_x_alt, y_a, 2, "qb_3"
        )
        qb4_refpoints = self._produce_qubit(
            float(self.coupler_widths[0]), qubits_center_x + qubit_spacing_x, y_b, 0, "qb_4"
        )
        qb5_refpoints = self._produce_qubit(
            float(self.coupler_widths[1]), qubits_center_x + qubit_spacing_x + qubit_spacing_x_alt, y_a, 2, "qb_5"
        )

        self._qubit_x_coords = [
            qubits_center_x - qubit_spacing_x,
            qubits_center_x - qubit_spacing_x + qubit_spacing_x_alt,
            qubits_center_x,
            qubits_center_x + qubit_spacing_x_alt,
            qubits_center_x + qubit_spacing_x,
            qubits_center_x + qubit_spacing_x_alt + qubit_spacing_x
        ]
        self._qubit_refpoints = [
            qb0_refpoints,
            qb1_refpoints,
            qb2_refpoints,
            qb3_refpoints,
            qb4_refpoints,
            qb5_refpoints
        ]

    def _produce_readout_resonator(self, capacitor, capacitor_dtrans, res_idx):

        total_length = float(self.readout_res_lengths[res_idx])
        turn_radius = 50

        # non-meandering part of the resonator
        meander_start = self.get_refpoints(capacitor, capacitor_dtrans)["port_a"]
        meander_end = self._qubit_refpoints[res_idx]["port_0"]

        # meandering part of the resonator
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

    def _produce_feedline_resonator(self):
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
        self._readout_structure_info["tees"] = [self.a, self.a, 2*self.a]
        self._tee_refpoints = tee_refpoints
####################
        feedline_tee_trans = pya.DTrans(
            0, False, 8600, 5e3
        )
        feedline_tee, feedline_tee_refp = self.insert_cell(cell_cross, feedline_tee_trans)
####################
        if isinstance(self.feedline_capacitor_n_fingers, tuple):
            self.feedline_capacitor_n_fingers = self.feedline_capacitor_n_fingers[0]
        cplr_params = cap_params(
            self.feedline_capacitor_n_fingers, float(self.l_fingers[0]), "smooth", finger_gap = self.b
        )
        cplr = self.add_element(**cplr_params)
        cplr_refpoints_rel = self.get_refpoints(cplr)
        cplr_dtrans = pya.DTrans(0, False, self.refpoints_fl_in["base"].x + 200, self.refpoints_fl_in["base"].y)
        _, cplr_refpoints = self.insert_cell(
            cplr, cplr_dtrans, rec_levels=None
        )

        self._produce_waveguide(
            [
                self.refpoints_fl_in["base"],
                cplr_refpoints["port_a"],
                # tee_refpoints[0]["port_left"],
            ]
        )
        self._produce_waveguide(
            [
                cplr_refpoints["port_b"],
                pya.DPoint(cplr_refpoints["port_b"].x + 200, cplr_refpoints["port_b"].y),
                # pya.DPoint(cplr_refpoints["port_b"].x + 80, 4530),
                pya.DPoint(cplr_refpoints["port_b"].x + 200, 5470),
                pya.DPoint(cplr_refpoints["port_b"].x + 200, 7500),
                pya.DPoint(cplr_refpoints["port_b"].x + 400, 7500),
                pya.DPoint(cplr_refpoints["port_b"].x + 400, 5470),
                pya.DPoint(cplr_refpoints["port_b"].x + 600, 5470),
                pya.DPoint(cplr_refpoints["port_b"].x + 600, 7500),
                pya.DPoint(cplr_refpoints["port_b"].x + 800, 7500),
                pya.DPoint(cplr_refpoints["port_b"].x + 800, 5e3),
                tee_refpoints[0]["port_left"],
            ], turn_radius= 50
        )
        
        self._produce_waveguide(
            [
                tee_refpoints[0]["port_right"],
                tee_refpoints[1]["port_right"],
            ]
        )
        
        self._produce_waveguide(
            [
                tee_refpoints[1]["port_left"],
                tee_refpoints[2]["port_left"],
            ]
        )
        self._produce_waveguide(
            [
                tee_refpoints[2]["port_right"],
                tee_refpoints[3]["port_right"],
            ]
        )
        self._produce_waveguide(
            [
                tee_refpoints[3]["port_right"],
                tee_refpoints[4]["port_right"],
            ]
        )
        self._produce_waveguide(
            [
                tee_refpoints[4]["port_right"],
                tee_refpoints[5]["port_right"],
            ]
        )
        self._produce_waveguide(
            [
                tee_refpoints[5]["port_left"],
                feedline_tee_refp["port_left"],
            ]
        )
        self._produce_waveguide(
            [
                feedline_tee_refp["port_right"],
                self.refpoints_fl_out["base"]
            ]
        )
        self._produce_waveguide(
            [
                feedline_tee_refp["port_bottom"],
                pya.DPoint(feedline_tee_refp["port_bottom"].x, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 100, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 200, feedline_tee_refp["port_bottom"].y - 500),
                pya.DPoint(feedline_tee_refp["port_bottom"].x + 200, feedline_tee_refp["port_bottom"].y - 1200),
            ]
        )

    def _produce_readout_resonators(self):
        # Coupler
        for i in range(6):
            cplr_params = cap_params(
                float(self.n_fingers[i]), float(self.l_fingers[i]), self.type_coupler[i], finger_gap = self.b
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
            
    # def get_readout_structure_info(self):
    #     if self.with_feedline_resonator:
    #         self._readout_structure_info["readout_structure"] = ["feedline_resonator"]
    #     else:
    #         self._readout_structure_info["readout_structure"] = ["simple"]

    #     df = pd.DataFrame.from_dict(self._readout_structure_info, orient = "index")
    #     df.to_csv("C:/Users/labranca/Desktop/work/QDAST/src/modeling/single_clockmons/single_clockmons_readout_structure.csv", mode='w+')
