# Import any modules, classes or functions used in our code.
from kqcircuits.elements.element import Element
from kqcircuits.pya_resolver import pya
from kqcircuits.util.parameters import Param, pdt
from kqcircuits.util.symmetric_polygons import polygon_with_vsym


# Any KQCircuits element must inherit from Element.
class SimpleCross(Element):

    # Define parameters for this class here.
    # Each parameter definition contains the parameter type, description and default value.
    # Other optional data such as the unit can also be defined for parameters.
    arm_length = Param(pdt.TypeDouble, "Cross arm length", 100, unit="μm")

    # The build() function is where the element geometry is built.
    def build(self):
        # We define a hardcoded value for arm_width, so it cannot be changed from outside like arm_length.
        arm_width = 30
        # Define some variables to hold values used commonly in this function.
        len1 = arm_width / 2
        len2 = arm_width / 2 + self.arm_length
        # Define the cross polygon using a list of DPoints.
        cross_poly = pya.DPolygon(
            [
                pya.DPoint(-len1, -len2),
                pya.DPoint(-len1, -len1),
                pya.DPoint(-len2, -len1),
                pya.DPoint(-len2, len1),
                pya.DPoint(-len1, len1),
                pya.DPoint(-len1, len2),
                pya.DPoint(len1, len2),
                pya.DPoint(len1, len1),
                pya.DPoint(len2, len1),
                pya.DPoint(len2, -len1),
                pya.DPoint(len1, -len1),
                pya.DPoint(len1, -len2),
            ]
        )
        # Add the cross polygon to the cell.
        # We use the get_layer() function to select in which layer the polygon is added.
        self.cell.shapes(self.get_layer("base_metal_gap_wo_grid")).insert(cross_poly)
