import os
from pathlib import Path

from kqcircuits import defaults
from kqcircuits.pya_resolver import pya
from kqcircuits.util.import_helper import module_from_file

SRC_PATH = Path(os.path.dirname(os.path.realpath(__file__)))
# workaround for Windows because os.path.realpath doesn't work there before Python 3.8
if os.name == "nt" and os.path.islink(Path(__file__).parent):
    SRC_PATH = Path(os.readlink(str(Path(__file__).parent)))
defaults.SRC_PATHS.append(SRC_PATH)
 
KQC_ROOT_PATH = Path(os.getenv("KQC_ROOT_PATH", defaults.ROOT_PATH))
defaults.ROOT_PATH = SRC_PATH.parent
defaults.TMP_PATH = defaults.TMP_PATH / "qdast" if os.getenv("KQC_TMP_PATH") else defaults.ROOT_PATH.joinpath("tmp")
defaults.TMP_PATH.mkdir(exist_ok=True)
defaults.SCRIPTS_PATH = defaults.ROOT_PATH.joinpath("scripts")
defaults.DRC_PATH = defaults.ROOT_PATH.joinpath("drc")
defaults.ANSYS_SCRIPT_PATHS += [defaults.SCRIPTS_PATH.joinpath("simulations").joinpath("ansys")]
defaults.ELMER_SCRIPT_PATHS += [defaults.SCRIPTS_PATH.joinpath("simulations").joinpath("elmer")]
defaults.VERSION_PATHS["QDAST"] = defaults.ROOT_PATH
defaults.default_drc_runset = "qdast_drc.lydrc"
