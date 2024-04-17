"""
Utilities to build interpolation libraries for Clockmon coupler capacitance data.

This module loads a CSV swept dataset (produced by Q3D simulations) and
provides several interpolating functions useful for circuit parametrization:

- clockmon_library: maps coupler width -> 3x3 capacitance matrix (deembedded)
- clockmon_coupling_libraries: provides interpolants for c_qr and c_sigma
  relations useful for design lookups
- clockmon_cqr_to_ground: maps c_qr -> C11 (coupler to ground) after deembed
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy import interpolate

from modeling.waveguides.libraries.waveguide_library import waveguide_library


_PATTERN = re.compile(r"[-+]?\d*\.\d+|\d+")
_DEFAULT_FILENAME = "sweeps/clockmon_capacitance_library_sim_q3d_results.csv"


def _load_csv_and_parse(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load the CSV file and return (coupler_widths, CMatrix).

    Parameters
    ----------
    file_path
        Path to the CSV file containing the sweep. Expected columns: first
        two columns (meta) and remaining 3x3 blocks per row.

    Returns
    -------
    coupler_widths : np.ndarray
        1D array of parsed numeric coupler widths.
    CMatrix : np.ndarray
        Array of shape (N, 3, 3) with the capacitance matrices for each sweep
        point.
    """
    df = pd.read_csv(file_path)

    # Parse the textual width column using the compiled regex
    coupler_widths_str = df["coupler_extent"].values
    coupler_widths = []
    for s in coupler_widths_str:
        m = _PATTERN.search(str(s))
        coupler_widths.append(float(m.group()) if m else np.nan)
    coupler_widths = np.array(coupler_widths, dtype=float)

    # Remaining columns form a flattened 3x3 matrix per row
    vals = df.iloc[:, 2:].values
    CMatrix = vals.reshape(-1, 3, 3)

    return coupler_widths, CMatrix


def clockmon_library(deembed: int = 200):
    """Return an interpolant mapping coupler width -> deembedded capacitance matrix.

    The CSV file is loaded from the same directory as this module under the
    default sweeps path. The returned object is a SciPy interp1d instance with
    axis=0 so calling ``lib(width)`` yields a (3,3) matrix.

    Parameters
    ----------
    deembed
        Deembedding distance in micrometers (used to subtract a waveguide
        contribution from C11).

    Returns
    -------
    library, coupler_widths, CMatrix
        - library: scipy.interpolate.interp1d mapping width -> (3,3) matrix
        - coupler_widths: ndarray of sweep widths
        - CMatrix: ndarray shape (N,3,3) of the deembedded capacitance matrices
    """
    base = Path(__file__).parent
    file_path = base / _DEFAULT_FILENAME

    coupler_widths, CMatrix = _load_csv_and_parse(str(file_path))

    wg_lib = waveguide_library()
    CMatrix = CMatrix.copy()
    CMatrix[:, 0, 0] -= wg_lib(deembed)

    library = interpolate.interp1d(coupler_widths, CMatrix, axis=0)
    return library, coupler_widths, CMatrix


def clockmon_coupling_libraries(deembed: int = 200):
    """Build interpolants relating coupler width, c_sigma and c_qr.

    Returns two interp1d objects:
    - coupler_width_given_c_qr: c_qr -> coupler_width
    - c_sigma_given_coupler_width: coupler_width -> c_sigma
    """
    base = Path(__file__).parent
    file_path = base / _DEFAULT_FILENAME

    coupler_widths, CMatrix = _load_csv_and_parse(str(file_path))

    CMatrix = CMatrix.copy()
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] -= wg_lib(deembed)

    c_sigmas, c_qrs = get_csigma_cqr(CMatrix)

    coupler_width_given_c_qr = interpolate.interp1d(c_qrs, coupler_widths)
    c_sigma_given_coupler_width = interpolate.interp1d(coupler_widths, c_sigmas)
    return coupler_width_given_c_qr, c_sigma_given_coupler_width


def clockmon_cqr_to_ground(deembed: int = 200):
    """Return an interpolant mapping c_qr -> C11 (coupler-to-ground) after deembed.

    Useful when selecting a coupler geometry from a target c_qr value.
    """
    base = Path(__file__).parent
    file_path = base / _DEFAULT_FILENAME

    coupler_widths, CMatrix = _load_csv_and_parse(str(file_path))
    CMatrix = CMatrix.copy()
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] -= wg_lib(deembed)

    _, c_qr = get_csigma_cqr(CMatrix)
    library = interpolate.interp1d(c_qr, CMatrix[:, 0, 0], axis=0)
    return library


def get_csigma_cqr(cmatrix: np.ndarray):
    """Compute c_sigma and c_qr from capacitance matrices.

    The input ``cmatrix`` may be a single 3x3 matrix or an array of shape
    (N,3,3). The function returns two arrays (or scalars) with the same
    leading dimension as the input.

    Returns
    -------
    c_sigma, c_qr
    """
    # Extract components (works for both 2D and 3D arrays)
    if cmatrix.ndim == 2:
        C12 = cmatrix[0, 1]
        C13 = cmatrix[0, 2]
        C22 = cmatrix[1, 1]
        C23 = cmatrix[1, 2]
        C33 = cmatrix[2, 2]
    else:
        C12 = cmatrix[:, 0, 1]
        C13 = cmatrix[:, 0, 2]
        C22 = cmatrix[:, 1, 1]
        C23 = cmatrix[:, 1, 2]
        C33 = cmatrix[:, 2, 2]

    # Formulas adapted from https://www.nature.com/articles/s41534-020-0269-1.pdf
    c_sigma = ((C33 + C13) * (C22 + C12)) / (C33 + C22 + C13 + C12) + C23
    beta = (C33 * C12 - C22 * C13) / (
        (C33 + C13) * (C22 + C12) + (C33 + C22 + C13 + C12) * C23
    )
    c_qr = beta * c_sigma
    return c_sigma, c_qr
