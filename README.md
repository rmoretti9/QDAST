# QubitDM

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)]()
[![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue)]()
[![Docs](https://img.shields.io/badge/docs-in_progress-lightgrey)]()

**Qubits for Direct Dark Matter Search**

This repository accompanies research from Roberto Moretti's PhD work. It provides tools and scripts for simulating interactions between dark matter candidates and superconducting transmon qubits, and for retrieving and analysing the outputs.

---

## Table of contents

- [Quick start](#quick-start)
- [Features](#features)
- [Installation](#installation)
- [Usage / Examples](#usage--examples)
- [Repository structure](#repository-structure)
- [Development](#development)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## Quick start

Clone the repository and install in editable mode:

```bash
git clone https://github.com/rmoretti9/QubitDM.git
cd QubitDM
pip install -e .
```

(If you use `conda` / virtualenv, create and activate your environment first.)

---

## Features

- Simulation modules for dark-matter ⇄ qubit interaction modeling
- Tools to run suites of simulations (runcards)
- Data analysis, management, and plotting helpers
- Example experiments and small demo scripts in `/examples`
---

## Installation

1. Create a virtual environment (recommended):

```bash
python -m venv .venv
source .venv/bin/activate   # Linux / macOS
.venv\Scripts\activate     # Windows (PowerShell)
```

2. Install the package and dependencies:

```bash
pip install -e .
# or
pip install -r requirements.txt
```

## Usage / Examples
1. Create a noise model
Configure a runcard in the ```runcards``` folder, then from ```QDAST```, run:
```bash
python src/qubitdm/create_noise_models.py --runcard_path <runcard_file_name>
```
Pass the argument ```--force``` to overwrite the noise model folder if existing.

This will create a folder in ```noise_models/V2/<filename>```. ```<filename>`` corresponds to the filename field you set in the runcard.

2. Simulate the dark matter search experiment with the noise model you just created
Always from ```QDAST```, run:
```bash
python src/qubitdm/simulate.py --sweep_folder noise_models/V2/<filename>
```
The results of the simulations will appear in ```noise_models/V2/<filename>/data_for_analysis.npy```, which can be analyzed with ```examples/sensitivity_study.ipynb```

3. Example
Run ```simulation.bat``` The ```.bat``` files help concatenating multiple experiments with different noise and experiment configurations.

## Repository structure

```
QubitDM/
├─ examples/                # runnable examples, notebooks and short tutorials
├─ src/                     # python package source (importable modules)
├─ runcards/                # runcards with the experimental settings and noise model data
├─ plots/                   # plots relevant to the analyses
├─ workload_results/        # workload .json files containing the outcome of circuit executions on real quantum hardware
├─ README.md
├─ requirements.txt
├─ pyproject.toml
└─ LICENSE
```

## Citation


---

## License

This project is released under the MIT License. See `LICENSE` for details.