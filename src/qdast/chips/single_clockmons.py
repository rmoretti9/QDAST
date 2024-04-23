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
class SingleClockmons(QDASTChip):
    """The PCell declaration for a SingleClockmons chip.

    The SingleClockmons chip has 4 qubits, which are coupled by lambda/2 readout resonators to the same feedline. The feedline
    crosses the center of the chip horizontally.  Half of the qubits are above the feedline and half are below it.

    Attributes:
        launchers: A dictionary where the keys are names of the launchers and values are tuples whose first elements
            are positions of the launchers.

        qubits_refpoints: A tuple where each element contains the refpoints for one of the qubits. The qubits are
            ordered such that 1,3 are the upper qubits (from left to right), while 0,2 are the lower qubits (from
            left to right).

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
    l_fingers = Param(
        pdt.TypeList,
        "Length of fingers for readout resonator couplers",
        [30, 30, 30, 30],
    )
    type_coupler = Param(
        pdt.TypeList,
        "Coupler type for test resonator couplers",
        ["smooth", "smooth", "smooth", "smooth"],
    )

    def build(self):
        """Produces a Clockmons PCell."""

        self.launchers = self.produce_launchers(
            "SquareNSWE_5x5", launcher_assignments={4: "FL-IN", 2: "FL-OUT"}
        )
        self.qubits_refpoints = self._produce_qubits()
        feedline_x_distance = 200
        self._produce_feedline(feedline_x_distance)
        self._produce_readout_resonators()

    def _produce_waveguide(self, path, term2=0, turn_radius=None):
        """Produces a coplanar waveguide that follows the given path.

        Args:
            path: a DPath object determining the waveguide path
            term2: term2 of the waveguide
            turn_radius: turn_radius of the waveguide

        Returns:
            length of the produced waveguide

        """
        if turn_radius is None:
            turn_radius = self.r
        waveguide = self.add_element(
            WaveguideCoplanar,
            path=pya.DPath(path, 1),
            r=turn_radius,
            term2=term2,
        )
        self.insert_cell(waveguide)
        return waveguide

    def _produce_qubit(self, qubit_cell, center_x, center_y, rotation, name=None):
        """Produces a qubit in a SingleClockmons chip.

        Args:
            qubit_cell: PCell of the qubit.
            center_x: X-coordinate of the center of the qubit.
            center_y: Y-coordinate of the center of the qubit.
            rotation: An integer which defines the rotation of the qubit in units of 90 degrees.
            name: A string containing the name of this qubit. Used to set the "id" property of the qubit instance.

        Returns:
            refpoints of the qubit.

        """
        qubit_trans = pya.DTrans(rotation, False, center_x, center_y)
        _, refpoints_abs = self.insert_cell(
            qubit_cell, qubit_trans, name, rec_levels=None
        )
        return refpoints_abs

    def _produce_qubits(self):
        """Produces four Clockmon qubits in predefined positions in a SingleClockmons chip.

        Three qubits are above the feedline and three are below it. The produced qubits are at equal distances
        from the feedline, and the distances between qubits in the x-direction are equal. The qubits are centered around
        the center of the chip.

        Returns:
            A tuple containing the refpoints of the qubits. Each element in the tuple contains the refpoints for a
            single qubit. The qubits are ordered such that 0,1 are the upper qubits (from left to right), while
            2,3 are the lower qubits (from left to right).

        """
        qubit = self.add_element(
            Clockmon,
            ground_gap=[650, 450],
            a=10,
            b=6,
            island_extent=[550, 130],
            coupler_extent=[150, 20],
            island_to_island_distance=20,
            coupler_offset=160,
            clock_diameter=60,
            bending_angle=0,
            junction_type="Manhattan Single Junction Centered",
        )
        qubit_spacing_x = 800  # shortest x-distance between qubit centers on different sides of the feedline
        qubit_spacing_y = 1500  # shortest y-distance between qubit centers on different sides of the feedline
        qubits_center_x = 2.5e3  # the x-coordinate around which qubits are centered
        # qubits above the feedline, from left to right
        y_a = 3.5e3 + qubit_spacing_y / 2
        y_b = 1.5e3 - qubit_spacing_y / 2
        qb0_refpoints = self._produce_qubit(
            qubit, qubits_center_x - qubit_spacing_x * (3 / 2), y_b, 0, "qb_0"
        )
        qb1_refpoints = self._produce_qubit(
            qubit, qubits_center_x - qubit_spacing_x * (1 / 2), y_a, 2, "qb_1"
        )
        # qubits below the feedline, from left to right
        qb2_refpoints = self._produce_qubit(
            qubit, qubits_center_x + qubit_spacing_x * (1 / 2), y_b, 0, "qb_2"
        )
        qb3_refpoints = self._produce_qubit(
            qubit, qubits_center_x + qubit_spacing_x * (3 / 2), y_a, 2, "qb_3"
        )
        self._qubit_x_coords = [
            qubits_center_x - qubit_spacing_x * (5 / 2),
            qubits_center_x - qubit_spacing_x * (3 / 2),
            qubits_center_x + qubit_spacing_x * (1 / 2),
            qubits_center_x + qubit_spacing_x * (3 / 2),
        ]
        self._qubit_refpoints = [
            qb0_refpoints,
            qb1_refpoints,
            qb2_refpoints,
            qb3_refpoints,
        ]
        return qb0_refpoints, qb1_refpoints, qb2_refpoints, qb3_refpoints

    def _produce_readout_resonator(self, capacitor, capacitor_dtrans, res_idx):

        total_length = float(self.readout_res_lengths[res_idx])
        turn_radius = 50

        # non-meandering part of the resonator
        meander_start = self.get_refpoints(capacitor, capacitor_dtrans)["port_a"]
        meander_end = self._qubit_refpoints[res_idx]["port_0"]

        # meandering part of the resonator
        w = 600
        num_meanders = _get_num_meanders(total_length, turn_radius, w)
        self.insert_cell(
            Meander,
            start=meander_start,
            end=meander_end,
            length=total_length,
            meanders=num_meanders,
            r=turn_radius,
        )

    def _produce_feedline(self, x_distance):
        """Produces a feedline for a SingleClockmons chip.

        The feedline is a waveguide connecting launcher "FL-IN" to launcher "FL-OUT". It goes horizontally from the
        launchers towards the center of the chip until it is clear of the junction test pads. Then it goes vertically
        to the center horizontal line, after which it follows the horizontal line until the two sides are connected.

        Args:
            x_distance: A float defining the x-distance of the vertical parts from the launchers. This is used to stay
                clear from the junction test pads.

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
                tee_refpoints[3]["port_left"],
                # pya.DPoint(self.launchers["FL-OUT"][0].x - x_distance, 2.5e3),
                pya.DPoint(
                    self.launchers["FL-OUT"][0].x - x_distance,
                    self.launchers["FL-OUT"][0].y,
                ),
                self.launchers["FL-OUT"][0],
            ]
        )

    def _produce_readout_resonators(self):
        # Coupler
        for i in range(4):
            cplr_params = cap_params(
                float(self.n_fingers[i]), float(self.l_fingers[i]), self.type_coupler[i]
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
