# QDAST — Quantum Design And Simulation Tool

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)

QDAST is a KQCircuits extension with custom PCells for planar superconducting qubit layouts and design automation.

> ⚠️ Compatibility: QDAST requires the forked KQCircuits repository maintained at `https://github.com/rmoretti9/KQCircuits.git`. Using other KQCircuits releases may break functionality. See **Requirements** for the exact commit/branch tested. :contentReference[oaicite:11]{index=11}

## Quick start

1. Install the required KQCircuits fork:
```bash
git clone https://github.com/rmoretti9/KQCircuits.git
python -m pip install -e klayout_package/python
python setup_within_klayout.py
```
2. Clone QDAST and install:
```bash
git clone https://github.com/rmoretti9/QDAST.git
cd QDAST
python -m pip install -e .
python setup_within_klayout.py
```
Then open KLayout → Libraries and verify that QDAST PCells appear.

## Project structure
```
QDAST/
├─ src/                       # package implementation (PCells & helpers)
├─ setup_within_klayout.py    # registration helper for KLayout
├─ requirements.txt
├─ pyproject.toml
├─ README.md
└─ LICENSE
```

## Cite

## Known issues
The package was originally developed for personal use and is not actively maintained. Use with caution; report issues via GitHub Issues if appropriate.

## License
This project is released under the MIT License. See `LICENSE` for details.