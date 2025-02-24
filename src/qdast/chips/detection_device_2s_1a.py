from math import pi, sqrt
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
from kqcircuits.elements.waveguide_coplanar_taper import WaveguideCoplanarTaper
from kqcircuits.util.geometry_helper import arc_points

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
@add_parameters_from(
    Clockmon,
    with_squid=False,
    sim_tool="eig"
)
class DetectionDevice2s1a(QDASTChip):
    readout_res_lengths = Param(
        pdt.TypeList,
        "Readout resonator lengths",
        [10000, 8100, 8200],
    )
    n_fingers = Param(
        pdt.TypeList,
        "Number of fingers for readout resonator couplers",
        [1.2, 1.3, 1.4],
    )
    feedline_capacitor_n_fingers = Param(
        pdt.TypeList,
        "Number of fingers for feedline couplers",
        [4, 4.3, 4.4],
    )
    def build(self):
        self.launchers = self.produce_launchers(
            "12-ports-10x10", launcher_assignments={}
        )
        self._cplr_res_refpoints = []
        self._produce_feedline_resonator(mirror_x = 1, mirror_y = 1, id = 0, l1_id = "9", l2_id = "10")
        self._produce_feedline_resonator(mirror_x = -1, mirror_y = 1, id = 1, l1_id = "7", l2_id = "6")
        self._produce_feedline_resonator(mirror_x = 1, mirror_y = -1, id = 2, l1_id = "1", l2_id = "12")

        self._produce_qubits()
        self._produce_readout_reasonators()

    def _produce_waveguide(self, path, term2=0, turn_radius=None, a = None, b = None, object = None):
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
            a = a,
            b = b
        )
        self.insert_cell(waveguide)
        if object == "feedline":
            self._readout_structure_info["feedline"].append(waveguide.length())
        return waveguide

    def _produce_qubits(self):
        # S1
        qubit = self.add_element(
            Clockmon, 
            ground_gap=[630, 610], 
            a=10,
            b=6,
            island_extent=[535, 200],
            coupler_widths=[0, 0, 0, 80, 0, 0],
            island_to_island_distance=50,
            coupler_offsets=[0, 0, 0, 255, 0, 0],
            clock_diameter=95,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered", 
            sim_tool = self.sim_tool,
            with_squid = self.with_squid,
            pad_width = 6,
            taper_width = 95/7
        )
        qubit_trans = pya.DCplxTrans(1, 0, False, 3500, 4000)
        _, refpoints_s1 = self.insert_cell(
            qubit, qubit_trans, name = "S1", rec_levels=None
        )
        # S2
        qubit = self.add_element(
            Clockmon, 
            ground_gap=[630, 610], 
            a=10,
            b=6,
            island_extent=[535, 200],
            coupler_widths=[0, 0, 0, 80, 0, 0],
            island_to_island_distance=50,
            coupler_offsets=[0, 0, 0, 255, 0, 0],
            clock_diameter=95,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered", 
            sim_tool = self.sim_tool,
            with_squid = self.with_squid,
            pad_width = 6,
            taper_width = 95/7
        )
        qubit_trans = pya.DCplxTrans(1, 0, False, 6500, 4000)
        _, refpoints_s2 = self.insert_cell(
            qubit, qubit_trans, name = "S2", rec_levels=None
        )
        # A1
        qubit = self.add_element(
            Clockmon, 
            ground_gap=[630, 610], 
            a=10,
            b=6,
            island_extent=[535, 200],
            coupler_widths=[80, 0, 0, 0, 0, 0],
            island_to_island_distance=50,
            coupler_offsets=[255, 0, 0, 0, 0, 0],
            clock_diameter=95,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered", 
            sim_tool = self.sim_tool,
            with_squid = self.with_squid,
            pad_width = 6,
            taper_width = 95/7
        )
        qubit_trans = pya.DCplxTrans(1, 0, False, 5000, 6000)
        _, refpoints_a1 = self.insert_cell(
            qubit, qubit_trans, name = "S1", rec_levels=None
        )
        self._qubits_refpoints = [refpoints_s1, refpoints_s2, refpoints_a1]


    def _produce_feedline_resonator(self, mirror_x, mirror_y, id, l1_id, l2_id):
        cell_cross = self.add_element(
            WaveguideCoplanarSplitter,
            **t_cross_parameters(
                a=self.a, b=self.b, a2=self.a, b2=self.b, length_extra_side=2 * self.a
            ),
        )
        rot_tee = 0 if mirror_y == 1 else 180

        feedline_tee_trans = pya.DCplxTrans(
            1, rot_tee, False, self.launchers[l2_id][0].x + mirror_x*250, self.launchers[l2_id][0].y
        )

        _, feedline_tee_refp = self.insert_cell(cell_cross, feedline_tee_trans)
        if mirror_x == 1 and mirror_y == 1:
            tee_port_1 = 'port_left'
            tee_port_2 = 'port_right'
        else:
            tee_port_1 = 'port_right'
            tee_port_2 = 'port_left'
        self._produce_waveguide(
            [
                self.launchers[l2_id][0],
                pya.DPoint(
                    feedline_tee_refp[tee_port_1].x,
                    self.launchers[l2_id][0].y,
                ),
            ]
        )
        self._produce_waveguide(
            [
                feedline_tee_refp[tee_port_2],
                pya.DPoint(
                    feedline_tee_refp[tee_port_2].x + mirror_x*50,
                    feedline_tee_refp[tee_port_2].y,
                ),
            ])
        self._produce_waveguide(
            [   
                pya.DPoint(
                    feedline_tee_refp[tee_port_2].x + mirror_x*750,
                    feedline_tee_refp[tee_port_2].y,
                ),
                pya.DPoint(
                    self.launchers[l1_id][0].x - mirror_x*500,
                    self.launchers[l2_id][0].y,
                ),
                pya.DPoint(
                    self.launchers[l1_id][0].x + mirror_x*(- 250 - (self.a/2 + self.b)*sqrt(2)/2),
                    self.launchers[l2_id][0].y + mirror_y*(- 250 + (self.a/2 + self.b)*sqrt(2)/2),
                )
            ]
        )
        if mirror_x == 1 and mirror_y == 1:
            rot_tee = 135
        if mirror_x == -1 and mirror_y == 1:
            rot_tee = 225
        if mirror_x == 1 and mirror_y == -1:
            rot_tee = 45
        resonator_tee_trans = pya.DCplxTrans(
            1, rot_tee, False, pya.DPoint(
                    self.launchers[l1_id][0].x - mirror_x*250,
                    self.launchers[l2_id][0].y - mirror_y*250,
                )
        )
        resonator_tee, resonator_tee_refp = self.insert_cell(cell_cross, resonator_tee_trans)

        # feedline capacitor
        cplr_params = cap_params(
            float(self.feedline_capacitor_n_fingers[id]), 0, "smooth", finger_gap = self.b
        )
        cplr = self.add_element(**cplr_params)
        if mirror_x == 1 and mirror_y == 1:
            rot_cplr = 90
        if mirror_x == -1 and mirror_y == 1:
            rot_cplr = 90
        if mirror_x == 1 and mirror_y == -1:
            rot_cplr = -90
        cplr_dtrans = pya.DCplxTrans(1, rot_cplr, False, self.launchers[l1_id][0].x, self.launchers[l1_id][0].y + mirror_y*300)
        _, cplr_refpoints = self.insert_cell(
            cplr, cplr_dtrans, rec_levels=None
        )

        # resonator capacitor
        cplr_res_params = cap_params(
            float(self.n_fingers[id]), 0, "smooth", finger_gap = self.b
        )
        cplr_res = self.add_element(**cplr_res_params)
        cplr_length  = self.get_refpoints(cplr_res)["port_b"].x

        if mirror_x == 1 and mirror_y == 1:
            rot_cplr = 45
        if mirror_x == -1 and mirror_y == 1:
            rot_cplr = 135
        if mirror_x == 1 and mirror_y == -1:
            rot_cplr = -45

        cplr_res_dtrans = pya.DCplxTrans(1, rot_cplr, False, pya.DPoint(resonator_tee_refp["port_bottom"].x + mirror_x*cplr_length*sqrt(2)/2, resonator_tee_refp["port_bottom"].y+ mirror_y*cplr_length*sqrt(2)/2))
        _, cplr_res_refpoints = self.insert_cell(
            cplr_res, cplr_res_dtrans, rec_levels=None
        )
        self._cplr_res_refpoints.append(cplr_res_refpoints)
        self._produce_waveguide(
            [
                resonator_tee_refp[tee_port_1],
                pya.DPoint(
                    self.launchers[l1_id][0].x,
                    self.launchers[l2_id][0].y- mirror_y*500,
                ),
                pya.DPoint(
                    self.launchers[l1_id][0].x,
                    cplr_refpoints['port_b'].y,
                ),
            ]
        )

        self._produce_waveguide([
                pya.DPoint(
                    self.launchers[l1_id][0].x,
                    cplr_refpoints['port_a'].y,
                ),
                self.launchers[l1_id][0]
        ])

        self._produce_waveguide([
                feedline_tee_refp['port_bottom'],
                pya.DPoint(
                    feedline_tee_refp['port_bottom'].x,
                    self.launchers[l2_id][0].y - mirror_y*500,
                ),
        ])

        # meander for tail
        meander_start = pya.DPoint(
                    feedline_tee_refp['port_bottom'].x,
                    self.launchers[l2_id][0].y - mirror_y*500,
                )
        meander_end = pya.DPoint(
                    feedline_tee_refp['port_bottom'].x,
                    self.launchers[l2_id][0].y - mirror_y*750,
                )
        w = 300
        meander_length = 500
        num_meanders = _get_num_meanders(meander_length, 50, w)
        self.insert_cell(
            Meander,
            start_point=meander_start,
            end_point=meander_end,
            length=meander_length,
            meanders=num_meanders,
            r=50,)
        

        # meander before split
        meander_start = pya.DPoint(
                    feedline_tee_refp[tee_port_2].x + mirror_x*50,
                    feedline_tee_refp[tee_port_2].y,
                )
        meander_end = pya.DPoint(
                    feedline_tee_refp[tee_port_2].x + mirror_x*750,
                    feedline_tee_refp[tee_port_2].y,
                )
        w = 300
        meander_length = 1200
        num_meanders = _get_num_meanders(meander_length, 50, w)
        self.insert_cell(
            Meander,
            start_point=meander_start,
            end_point=meander_end,
            length=meander_length,
            meanders=num_meanders,
            r=50,)
    def _produce_readout_reasonators(self):
        wg_0 = self._produce_waveguide([
            self._cplr_res_refpoints[0]['port_b'],
            pya.DPoint(
                self._cplr_res_refpoints[0]['port_b'].x + 200,
                self._cplr_res_refpoints[0]['port_b'].y + 200,
            ),
            pya.DPoint(
                self._cplr_res_refpoints[0]['port_b'].x + 300,
                self._cplr_res_refpoints[0]['port_b'].y + 300,
            ),
            pya.DPoint(
                self._cplr_res_refpoints[0]['port_b'].x + 800,
                self._cplr_res_refpoints[0]['port_b'].y - 200,
            ),
            pya.DPoint(
                self._qubits_refpoints[0]['port_3'].x,
                self._cplr_res_refpoints[0]['port_b'].y - 200,
            ),
            pya.DPoint(
                self._qubits_refpoints[0]['port_3'].x,
                self._cplr_res_refpoints[0]['port_b'].y - 100,
            ),
        ])
        meander_start = pya.DPoint(
                    self._qubits_refpoints[0]['port_3'].x,
                    self._cplr_res_refpoints[0]['port_b'].y -100,
                )
        meander_end = self._qubits_refpoints[0]['port_3']
        w = 1000
        meander_length = float(self.readout_res_lengths[0]) - wg_0.length()
        num_meanders = _get_num_meanders(meander_length, 50, w)
        self.insert_cell(
            Meander,
            start_point=meander_start,
            end_point=meander_end,
            length=meander_length,
            meanders=num_meanders,
            r=50,)
        


        wg_1 = self._produce_waveguide([
            self._cplr_res_refpoints[1]['port_b'],
            pya.DPoint(
                self._cplr_res_refpoints[1]['port_b'].x - 200,
                self._cplr_res_refpoints[1]['port_b'].y + 200,
            ),
            pya.DPoint(
                self._cplr_res_refpoints[1]['port_b'].x - 300,
                self._cplr_res_refpoints[1]['port_b'].y + 300,
            ),
            pya.DPoint(
                self._cplr_res_refpoints[1]['port_b'].x - 800,
                self._cplr_res_refpoints[0]['port_b'].y - 200,
            ),
            pya.DPoint(
                self._qubits_refpoints[1]['port_3'].x,
                self._cplr_res_refpoints[1]['port_b'].y - 200,
            ),
            pya.DPoint(
                self._qubits_refpoints[1]['port_3'].x,
                self._cplr_res_refpoints[1]['port_b'].y - 100,
            ),
        ])
        meander_start = pya.DPoint(
                    self._qubits_refpoints[1]['port_3'].x,
                    self._cplr_res_refpoints[1]['port_b'].y -100,
                )
        meander_end = self._qubits_refpoints[1]['port_3']
        w = 1000
        meander_length = float(self.readout_res_lengths[1]) - wg_1.length()
        num_meanders = _get_num_meanders(meander_length, 50, w)
        self.insert_cell(
            Meander,
            start_point=meander_start,
            end_point=meander_end,
            length=meander_length,
            meanders=num_meanders,
            r=50,)