import os
from pathlib import Path

from kqcircuits import defaults
from kqcircuits.pya_resolver import pya

SRC_PATH = Path(os.path.dirname(os.path.realpath(__file__)))
# workaround for Windows because os.path.realpath doesn't work there before Python 3.8
if os.name == "nt" and os.path.islink(Path(__file__).parent):
    SRC_PATH = Path(os.readlink(str(Path(__file__).parent)))
defaults.SRC_PATHS.append(SRC_PATH)

KQC_ROOT_PATH = Path(os.getenv("KQC_ROOT_PATH", defaults.ROOT_PATH))
defaults.ROOT_PATH = SRC_PATH.parent
defaults.TMP_PATH = (
    defaults.TMP_PATH / "qdast"
    if os.getenv("KQC_TMP_PATH")
    else defaults.ROOT_PATH.joinpath("tmp")
)
defaults.TMP_PATH.mkdir(exist_ok=True)
defaults.ANSYS_SCRIPT_PATHS += [
    defaults.SCRIPTS_PATH.joinpath("simulations").joinpath("ansys")
]
defaults.ELMER_SCRIPT_PATHS += [
    defaults.SCRIPTS_PATH.joinpath("simulations").joinpath("elmer")
]
defaults.VERSION_PATHS["QDAST"] = defaults.ROOT_PATH

default_sampleholders = {
    "SquareNSWE_5x5": {
        "n": 4,
        "launcher_type": "RF",
        "launcher_width": 240,
        "launcher_gap": 144,
        "launcher_indent": 680,
        "pad_pitch": 1200,
        "chip_box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(5000, 5000)),
    },

    "6-ports-10x10": {
        "n": 6,
        "launcher_type": "RF",
        "launcher_width": 200,
        "launcher_gap": 153,
        "launcher_indent": 680,
        "pad_pitch": 1200,
        "chip_box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(10000, 10000)),
    },


    "4-ports-75x75": {
        "n": 4,
        "launcher_type": "RF",
        "launcher_width": 200,
        "launcher_gap": 153,
        "launcher_indent": 773,
        "pad_pitch": 1200,
        "chip_box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(7500, 7500)),
    },

    "12-ports-10x10": {
        "n": 12,
        "launcher_type": "RF",
        "launcher_width": 250,
        "launcher_gap": 140,
        "launcher_indent": 800,
        "launcher_frame_gap": 180,
        "pad_pitch": 2500,
        "chip_box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(10000, 10000)),
    },
}