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
link_map = (
    ("src/qdast", "python/qdast"),)

setup_symlinks(qdast_root_path, configdir, link_map, args.unlink)

# # TODO install required packages: kqcircuits

# if (
#     not args.unlink and os.name == "posix" and sys.platform == "darwin"
# ):  # Install missing dependecy in MacOS
#     cmd = "pip install requests>=2.27.1"
#     td = get_klayout_packages_path(
#         "/Applications/klayout.app/Contents/Frameworks/Python.framework/Versions"
#     )
#     if not os.path.exists(td):
#         # Homebrew installs under /Applications/Klayout/klayout.app
#         td = get_klayout_packages_path(
#             "/Applications/KLayout/klayout.app/Contents/Frameworks/Python.framework/Versions"
#         )
#     # KLayout may use either its own site-packages or the system site-packages, depending on the build
#     if os.path.exists(td):
#         cmd += f' --target="{td}"'
#     os.system(cmd)

print("Finished setting up QDAST.")
