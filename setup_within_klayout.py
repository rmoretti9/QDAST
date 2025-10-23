import os
import sys
import argparse

# Check configured KQC at .klayout
if os.name == "nt":
    config_dir_name = "KLayout"
elif os.name == "posix":
    config_dir_name = ".klayout"
else:
    raise SystemError("Error: unsupported operating system")

parser = argparse.ArgumentParser(description="QDAST setup within klayout")
parser.add_argument(
    "alt_name",
    nargs="?",
    type=str,
    help="Alternative installation config directory name",
)
parser.add_argument("--unlink", action="store_true", help="remove links")
args = parser.parse_args()

configdir = os.path.join(os.path.expanduser("~"), config_dir_name)
if args.alt_name is not None:
    configdir = f"{configdir}_alt/{args.alt_name}"

kqc = os.path.join(configdir, "python", "kqcircuits")
if not os.path.exists(kqc):
    print("Please run setup_within_klayout.py in KQC first")
    sys.exit(-1)

sys.path.append(f"{os.path.realpath(kqc)}/../../..")  # KQC's root directory
from setup_helper import setup_symlinks

# QDAST source path
qdast_root_path = os.path.dirname(os.path.abspath(__file__))

# create symlink between KLayout python folder and qdast folder
link_map = (("src/qdast", "python/qdast"),)

setup_symlinks(qdast_root_path, configdir, link_map, args.unlink)
print("Finished setting up QDAST.")
