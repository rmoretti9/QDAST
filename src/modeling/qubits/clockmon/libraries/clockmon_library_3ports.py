"""
Utilities for 3-port Clockmon capacitance libraries.

This module loads Q3D CSV sweeps for 3-port clockmon geometries and builds
interpolators mapping coupler geometries (or derived quantities like c_qr)
onto capacitance-matrix entries or geometric parameters.

The refactor keeps the original numeric logic, factors CSV parsing, adds
clearer docstrings and light typing hints while preserving API and outputs.
"""
from __future__ import annotations

import ast
import os
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy import interpolate

from modeling.waveguides.libraries.waveguide_library import waveguide_library
from qucat import Network, L, J, C

_DEFAULT_SUFFIX = "_sim_q3d_results.csv"


def _load_5x5_csv(file_path: str) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Load a CSV produced for 3-port clockmon sweeps.

    Returns (parsed_cplr_lists, CMatrix, dataframe).
    CMatrix has shape (N,5,5).
    """
    df = pd.read_csv(file_path)
    parsed = df["coupler_widths"].apply(lambda x: ast.literal_eval(x)).values

    N = len(df)
    CMatrix = np.zeros((N, 5, 5))
    for i in range(N):
        CMatrix[i] = df.iloc[i].values[2:].reshape(5, 5)

    return parsed, CMatrix, df


def clockmon_library_3ports(deembed: int = 200, port_id: str = "3_4_5") -> interpolate.LinearNDInterpolator:
    """Return a 3D interpolator mapping three coupler widths -> 5x5 capacitance.

    Parameters
    ----------
    deembed
        Deembedding distance in micrometers (subtracted from three diagonal
        capacitances using the waveguide library).
    port_id
        String of the form "i_j_k" selecting which coupler indexes correspond
        to the three ports present in the CSV (default "3_4_5").

    Returns
    -------
    interpolator
        ``scipy.interpolate.LinearNDInterpolator`` that maps
        (cplr0_width, cplr1_width, cplr2_width) -> 5x5 capacitance matrix.
    """
    base = Path(__file__).parent
    file_name = f"sweeps/clockmon_capacitance_library_{port_id}{_DEFAULT_SUFFIX}"
    file_path = base / file_name

    parsed_lists, CMatrix, _ = _load_5x5_csv(str(file_path))

    i0, i1, i2 = [int(x) for x in port_id.split("_")]
    cplr0_widths = np.array([float(lst[i0]) for lst in parsed_lists])
    cplr1_widths = np.array([float(lst[i1]) for lst in parsed_lists])
    cplr2_widths = np.array([float(lst[i2]) for lst in parsed_lists])

    wg_lib = waveguide_library()
    CMatrix = CMatrix.copy()
    CMatrix[:, 0, 0] -= wg_lib(deembed)
    CMatrix[:, 1, 1] -= wg_lib(deembed)
    CMatrix[:, 2, 2] -= wg_lib(deembed)

    library = interpolate.LinearNDInterpolator((cplr0_widths, cplr1_widths, cplr2_widths), CMatrix)
    return library


def clockmon_coupling_libraries(deembed: int = 200, port_id: str = "3_4_5") -> interpolate.LinearNDInterpolator:
    """Return an interpolator mapping (c_qr_0,c_qr_1,c_qr_2) -> (w0,w1,w2).

    Useful to pick coupler geometries from target effective couplings.
    """
    base = Path(__file__).parent
    file_name = f"sweeps/clockmon_capacitance_library_{port_id}{_DEFAULT_SUFFIX}"
    file_path = base / file_name

    parsed_lists, CMatrix, _ = _load_5x5_csv(str(file_path))

    i0, i1, i2 = [int(x) for x in port_id.split("_")]
    cplr0_widths = np.array([float(lst[i0]) for lst in parsed_lists])
    cplr1_widths = np.array([float(lst[i1]) for lst in parsed_lists])
    cplr2_widths = np.array([float(lst[i2]) for lst in parsed_lists])

    CMatrix = CMatrix.copy()
    wg_lib = waveguide_library()
    CMatrix[:, 0, 0] -= wg_lib(deembed)
    CMatrix[:, 1, 1] -= wg_lib(deembed)
    CMatrix[:, 2, 2] -= wg_lib(deembed)

    c_qr_0, c_qr_1, c_qr_2 = get_cqr(CMatrix)

    points = np.column_stack((c_qr_0, c_qr_1, c_qr_2))
    values = np.column_stack((cplr0_widths, cplr1_widths, cplr2_widths))

    library = interpolate.LinearNDInterpolator(points, values)
    return library


def get_csigma_cqr(CMatrix: np.ndarray):
    """Compute c_sigma and the three c_qr values for a 5x5 capacitance matrix.

    The routine builds a lumped Network (qucat) with a Josephson inductance
    to extract a mode frequency and infer c_sigma from the mode.

    Returns
    -------
    c_sigma, c_qr_0, c_qr_1, c_qr_2
    """
    dim = CMatrix.shape[0]
    network = []
    for i in range(dim + 1):
        for j in range(i + 1, dim + 1):
            if i == 0:
                network.append(C(i, j, CMatrix[j - 1, j - 1]))
            else:
                network.append(C(i, j, CMatrix[i - 1, j - 1]))
    network.append(L(dim - 1, dim, "Lj"))
    cir = Network(network)
    Lj = 15e-9

    # use the qucat helper that returns eigenfrequencies
    f = cir.eigenfrequencies(Lj=Lj)
    c_sigma = 1 / Lj / (2 * np.pi * max(f)) ** 2
    c_qr_0, c_qr_1, c_qr_2 = get_cqr(CMatrix)
    return c_sigma, c_qr_0, c_qr_1, c_qr_2


def get_cqr(CMatrix: np.ndarray):
    """Compute three effective couplings c_qr_i from a 5x5 capacitance matrix.

    Works with both shape (5,5) and (N,5,5) inputs and returns arrays of
    matching leading dimension when applicable.
    """
    if CMatrix.ndim == 2:
        C14 = CMatrix[0, 3]
        C15 = CMatrix[0, 4]
        C24 = CMatrix[1, 3]
        C25 = CMatrix[1, 4]
        C34 = CMatrix[2, 3]
        C35 = CMatrix[2, 4]
        C44 = CMatrix[3, 3]
        C55 = CMatrix[4, 4]
    else:
        C14 = CMatrix[:, 0, 3]
        C15 = CMatrix[:, 0, 4]
        C24 = CMatrix[:, 1, 3]
        C25 = CMatrix[:, 1, 4]
        C34 = CMatrix[:, 2, 3]
        C35 = CMatrix[:, 2, 4]
        C44 = CMatrix[:, 3, 3]
        C55 = CMatrix[:, 4, 4]

    c_qr_1 = np.abs(C14 * C55 - C15 * C44) / (C44 + C14 + C55 + C15)
    c_qr_2 = np.abs(C24 * C55 - C25 * C44) / (C44 + C24 + C55 + C25)
    c_qr_3 = np.abs(C34 * C55 - C35 * C44) / (C44 + C34 + C55 + C35)
    return c_qr_1, c_qr_2, c_qr_3


def clockmon_cqr_to_ground(deembed: int = 200, port_id: str = "3_4_5", ground_id: int = 0) -> interpolate.LinearNDInterpolator:
    """Return an interpolator mapping (c_qr_0,c_qr_1,c_qr_2) -> C_gg (ground diag).

    ground_id selects which diagonal element of the 5x5 matrix to return.
    """
    base = Path(__file__).parent
    file_name = f"sweeps/clockmon_capacitance_library_{port_id}{_DEFAULT_SUFFIX}"
    file_path = base / file_name

    parsed_lists, CMatrix, _ = _load_5x5_csv(str(file_path))
    wg_lib = waveguide_library()
    CMatrix = CMatrix.copy()
    CMatrix[:, 0, 0] -= wg_lib(deembed)
    CMatrix[:, 1, 1] -= wg_lib(deembed)
    CMatrix[:, 2, 2] -= wg_lib(deembed)

    c_qr_0, c_qr_1, c_qr_2 = get_cqr(CMatrix)
    library = interpolate.LinearNDInterpolator((c_qr_0, c_qr_1, c_qr_2), CMatrix[:, ground_id, ground_id])
    return library