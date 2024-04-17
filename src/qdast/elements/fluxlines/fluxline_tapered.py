from kqcircuits.pya_resolver import pya
from kqcircuits.elements.fluxlines.fluxline import Fluxline


class FluxlineTapered(Fluxline):
    """Fluxline variant "tapered"."""

    def build(self):
        # shorthands
        a = self.a  # waveguide center width
        b = self.b  # waveguide gap width
        taper_length = 100
        a2 = 4
        b2 = 25
        right_gap = pya.DPolygon(
            [
                pya.DPoint(a / 2, -50 - taper_length),
                pya.DPoint(a / 2 + b, -50 - taper_length),
                pya.DPoint(a / 2 + b, -taper_length),
                pya.DPoint(a2 / 2 + b2, 0),
                pya.DPoint(a2 / 2, 0),
                pya.DPoint(a / 2, -taper_length),
            ]
        )
        left_gap = pya.DPolygon(
            [
                pya.DPoint(-a / 2, -50 - taper_length),
                pya.DPoint(-a / 2 - b, -50 - taper_length),
                pya.DPoint(-a / 2 - b, -taper_length),
                pya.DPoint(-a2 / 2 - b2, 0),
                pya.DPoint(-a2 / 2, 0),
                pya.DPoint(-a / 2, -taper_length),
            ]
        )

        self._insert_fluxline_shapes(left_gap, right_gap)
        self._add_fluxline_refpoints(pya.DPoint(0, -50 - taper_length))
