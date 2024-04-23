from kqcircuits.elements.chip_frame import ChipFrame
from kqcircuits.util.parameters import Param, pdt
from kqcircuits.util.label import produce_label
from kqcircuits.pya_resolver import pya
from kqcircuits.util.label import produce_label, LabelOrigin
from kqcircuits.util.parameters import Param, pdt
from kqcircuits.defaults import default_chip_label_face_prefixes


class QDASTChipFrame(ChipFrame):
    """Wrapper of KQCircuits ChipFrame to further customize labeling."""
    produce_labels = Param(pdt.TypeBoolean, "Produce chip labels", True)

    def build(self):
        """Produces dicing edge, markers, labels and ground grid for the chip face."""
        self._produce_dicing_edge()
        if self.produce_labels:
            self._produce_labels()
        self._produce_markers()

    def _produce_labels(self):
        x_min, x_max, y_min, y_max = self._box_points()
        if self.use_face_prefix:
            face_id = self.face()["id"]
            face_prefix = (
                default_chip_label_face_prefixes[face_id].upper()
                if default_chip_label_face_prefixes
                and (face_id in default_chip_label_face_prefixes)
                else face_id.upper()
            )
            chip_name = face_prefix + self.name_chip
        else:
            chip_name = self.name_chip
        labels = [self.name_mask, chip_name, self.name_copy, self.name_brand]
        if labels[0] != "":
            self._produce_label(
                labels[0], pya.DPoint(x_min, y_max), LabelOrigin.TOPLEFT
            )
        if labels[1] != "":
            self._produce_label(
                labels[1], pya.DPoint(x_max, y_max), LabelOrigin.TOPRIGHT
            )
        if labels[2] != "":
            self._produce_label(
                labels[2], pya.DPoint(x_max, y_min), LabelOrigin.BOTTOMRIGHT
            )
        if labels[3] != "":
            self._produce_label(
                labels[3], pya.DPoint(x_min, y_min), LabelOrigin.BOTTOMLEFT
            )

    def _produce_label(self, label, location, origin):
        """Produces Text PCells with text `label` with `origin` of the text at `location`.

        Wrapper for the stand alone function `produce_label`.
        Text size scales with chip dimension for chips smaller than 7 mm.

        Args:
            label: the produced text
            location: DPoint of the location of the text
            origin: LabelOrigin of the corner of the label to be placed at the location

        Effect:
            label PCells added to the layout into the parent PCell
        """
        size = 200 * min(
            1, self.box.width() / 7000, self.box.height() / 7000
        )  # we might reduce further for smaller sample holders
        produce_label(
            self.cell,
            label,
            location,
            origin,
            self.dice_width,
            self.text_margin,
            [
                self.face()["base_metal_gap_wo_grid"],
                self.face()["base_metal_gap_for_EBL"],
            ],
            self.face()["ground_grid_avoidance"],
            size,
        )
