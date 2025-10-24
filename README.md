# QDAST

**QDAST – Quantum Design and Simulation Tool**  
A KQCircuits-based extension for planar superconducting qubit designs, providing custom PCells for quantum device layouts.

---

## Overview

QDAST is designed as an extension to **KQCircuits**, offering additional **custom PCells** for planar superconducting qubit design and simulation.

> **Warning:** This package is compatible only with the **forked version of KQCircuits** available at `https://github.com/rmoretti9/KQCircuits.git`. Using newer or different versions of KQCircuits may break QDAST functionality.  
> Note: QDAST was originally developed for personal use and **will not be actively maintained**.

---

## Installation

1. **Install KQCircuits**  
   Follow the Developer Standalone module setup instructions provided by the KQCircuits developers at `https://iqm-finland.github.io/KQCircuits/developer/standalone.html`.

2. **Clone QDAST repository**  
   ```bash
   git clone <repository_link>
   cd QDAST
    ```

3. **Install QDAST**  
   ```bash
    python -m pip install -e .
    python setup_within_klayout.py
    ```

4. **Verify Installation**  
Open KLayout and check the libraries. You should see that QDAST's PCells have been added.
Reloading KQCircuits libraries will also reload QDAST’s libraries automatically.