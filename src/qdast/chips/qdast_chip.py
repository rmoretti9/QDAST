from qdast.defaults import default_sampleholders
from qdast.elements.qdast_chip_frame import QDASTChipFrame
from kqcircuits.chips.chip import Chip
from kqcircuits.defaults import default_marker_type, default_layers
from kqcircuits.util.parameters import add_parameters_from, add_parameter
from kqcircuits.pya_resolver import pya


@add_parameter(QDASTChipFrame, "box", hidden=True)
@add_parameters_from(
    QDASTChipFrame,
    "name_mask",
    "name_chip",
    "name_copy",
    "name_brand",
    "chip_dicing_in_base_metal",
    "dice_grid_margin",
    "produce_labels",
    marker_types=[default_marker_type] * 8,
)
class QDASTChip(Chip):
    """Wrapper of KQCircuits Chip for further customization."""
    def produce_launchers(
        self, sampleholder_type, launcher_assignments=None, enabled=None, face_id=0
    ):
        """Produces launchers for typical sample holders and sets chip size (``self.box``) accordingly.

        This is a wrapper around ``produce_n_launchers()`` to generate typical launcher configurations.

        Args:
            sampleholder_type: name of the sample holder type
            launcher_assignments: dictionary of (port_id: name) that assigns a name to some of the launchers
            enabled: list of enabled launchers, empty means all
            face_id: index of face_ids in which to insert the launchers

        Returns:
            launchers as a dictionary :code:`{name: (point, heading, distance from chip edge)}`

        """

        if (
            sampleholder_type == "SMA8"
        ):  # this is special: it has default launcher assignments
            if not launcher_assignments:
                launcher_assignments = {
                    1: "NW",
                    2: "NE",
                    3: "EN",
                    4: "ES",
                    5: "SE",
                    6: "SW",
                    7: "WS",
                    8: "WN",
                }

        if sampleholder_type == "6-ports-10x10":
            return {}

        if sampleholder_type in default_sampleholders:
            return self.produce_n_launchers(
                **default_sampleholders[sampleholder_type],
                launcher_assignments=launcher_assignments,
                enabled=enabled,
                face_id=face_id,
            )
        return {}

    def produce_frame(self, frame_parameters, trans=pya.DTrans()):
        """Produces a chip frame and markers for the given face.

        Args:
            frame_parameters: PCell parameters for the chip frame
            trans: DTrans for the chip frame, default=pya.DTrans()
        """
        self.insert_cell(QDASTChipFrame, trans, **frame_parameters)