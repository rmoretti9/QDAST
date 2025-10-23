"""
Utilities for 2-port Clockmon capacitance libraries.

This module provides helpers that load Q3D CSV sweeps for 2-port
clockmon geometries and build interpolators mapping geometrical parameters
(or derived quantities like c_qr) to capacitance matrix entries.
"""
from __future__ import annotations

import ast
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy import interpolate

from modeling.waveguides.libraries.waveguide_library import waveguide_library
from qucat import Network, J, C

_PATTERN = r"[-+]?\d*\.\d+|\d+"  # kept for possible future use
_DEFAULT_SUFFIX = "_sim_q3d_results.csv"


def _load_4x4_csv(file_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load a CSV produced for 2-port clockmon sweeps.

    Returns (cplr0_widths, cplr2_widths, CMatrix).
    ``CMatrix`` has shape (N,4,4).
    """
    df = pd.read_csv(file_path)

    # "coupler_widths" column contains a Python list literal per row; parse it
    parsed = df["coupler_widths"].apply(lambda x: ast.literal_eval(x))
    # first and second coupler widths (indexing depends on the port_id caller)
    # here the caller will select appropriate indices
    cplr_lists = parsed.values

    N = len(df)
    CMatrix = np.zeros((N, 4, 4))
    for i in range(N):
        CMatrix[i] = df.iloc[i].values[2:].reshape(4, 4)

    return cplr_lists, CMatrix, df


def clockmon_library_2ports(deembed: int = 200, port_id: str = "0_2", alternative_geometry: bool = False):
    """Return a 2D interpolator mapping two coupler widths -> 4x4 capacitance.

    Parameters
    ----------
    deembed
        Deembedding distance in micrometers (subtracted from the diagonal
        C11 and C22 entries using the waveguide library).
    port_id
        String of the form "i_j" selecting which coupler indexes correspond
        to the two ports present in the CSV (default "0_2").
    alternative_geometry
        Use the "_sensing" variant of the CSV filename when True.

    Returns
    -------
    interpolator
        A ``scipy.interpolate.LinearNDInterpolator`` that maps
        (cplr0_width, cplr2_width) -> 4x4 capacitance matrix.
    """
    base = Path(__file__).parent
    pid = port_id
    suffix = "_sensing" if alternative_geometry else ""
    file_name = f"sweeps/clockmon_capacitance_library_{pid}{suffix}{_DEFAULT_SUFFIX}"
    file_path = base / file_name

    cplr_lists, CMatrix, df = _load_4x4_csv(str(file_path))

    # Extract the two coupler widths per row according to port_id
    i0, i1 = [int(x) for x in port_id.split("_")]
    cplr0_widths = np.array([float(lst[i0]) for lst in cplr_lists])
    cplr2_widths = np.array([float(lst[i1]) for lst in cplr_lists])

    # Deembed the waveguide contribution from the two coupling nodes
    wg_lib = waveguide_library()
    CMatrix = CMatrix.copy()
    CMatrix[:, 0, 0] -= wg_lib(deembed)
    CMatrix[:, 1, 1] -= wg_lib(deembed)

    library = interpolate.LinearNDInterpolator((cplr0_widths, cplr2_widths), CMatrix)
    return library


def clockmon_coupling_libraries(deembed: int = 200, port_id: str = "0_2", alternative_geometry: bool = False):
    """Build an interpolator mapping (c_qr_1, c_qr_2) -> (cplr0_width, cplr2_width).

    This is useful to pick coupler geometries from target effective couplings.
    """
    base = Path(__file__).parent
    pid = port_id
    suffix = "_sensing" if alternative_geometry else ""
    file_name = f"sweeps/clockmon_capacitance_library_{pid}{suffix}{_DEFAULT_SUFFIX}"
    file_path = base / file_name

    cplr_lists, CMatrix, df = _load_4x4_csv(str(file_path))

    # select indices (first entry is always parsed[0])
    cplr0_widths = np.array([float(lst[0]) for lst in cplr_lists])
    idx = int(port_id.split("_")[-1])
    cplr2_widths = np.array([float(lst[idx]) for lst in cplr_lists])

    # Deembed
    wg_lib = waveguide_library()
    CMatrix = CMatrix.copy()
    CMatrix[:, 0, 0] -= wg_lib(deembed)
    CMatrix[:, 1, 1] -= wg_lib(deembed)

    c_qr_1, c_qr_2 = get_cqr(CMatrix)

    points = np.column_stack((c_qr_1, c_qr_2))
    values = np.column_stack((cplr0_widths, cplr2_widths))

    library = interpolate.LinearNDInterpolator(points, values)
    return library


def get_csigma_cqr(CMatrix: np.ndarray):
    """Compute c_sigma and the two c_qr values for a 4x4 capacitance matrix.

    The function returns (c_sigma, c_qr_1, c_qr_2).
    """
    # Build a lumped network and extract the mode frequency to infer c_sigma
    dim = CMatrix.shape[0]
    network = []
    for i in range(dim + 1):
        for j in range(i + 1, dim + 1):
            if i == 0:
                network.append(C(i, j, CMatrix[j - 1, j - 1]))
            else:
                network.append(C(i, j, CMatrix[i - 1, j - 1]))
    # add a Josephson inductance placeholder
    network.append(J(3, 4, "Lj"))
    cir = Network(network)
    Lj = 10e-9
    f, _, _, _ = cir.f_k_A_chi(Lj=Lj)
    c_sigma = 1 / Lj / (2 * np.pi * f[0]) ** 2
    c_qr_1, c_qr_2 = get_cqr(CMatrix)
    return c_sigma, c_qr_1, c_qr_2


def get_cqr(CMatrix: np.ndarray):
    """Compute the two effective couplings c_qr_1 and c_qr_2 from a 4x4 matrix.

    Works with both shape (4,4) and (N,4,4) inputs and returns arrays of
    matching leading dimension when applicable.
    """
    if CMatrix.ndim == 2:
        C12 = CMatrix[0, 1]
        C13 = CMatrix[0, 2]
        C14 = CMatrix[0, 3]
        C22 = CMatrix[1, 1]
        C23 = CMatrix[1, 2]
        C24 = CMatrix[1, 3]
        C33 = CMatrix[2, 2]
        C34 = CMatrix[2, 3]
        C44 = CMatrix[3, 3]
    else:
        C12 = CMatrix[:, 0, 1]
        C13 = CMatrix[:, 0, 2]
        C14 = CMatrix[:, 0, 3]
        C22 = CMatrix[:, 1, 1]
        C23 = CMatrix[:, 1, 2]
        C24 = CMatrix[:, 1, 3]
        C33 = CMatrix[:, 2, 2]
        C34 = CMatrix[:, 2, 3]
        C44 = CMatrix[:, 3, 3]

    c_qr_1 = np.abs(C13 * C44 - C14 * C33) / (C33 + C13 + C44 + C14)
    c_qr_2 = np.abs(C23 * C44 - C24 * C33) / (C33 + C23 + C44 + C24)
    return c_qr_1, c_qr_2


def clockmon_cqr_to_ground(deembed: int = 200, port_id: str = "0_2", ground_id: int = 0, alternative_geometry: bool = False):
    """Return an interpolator mapping (c_qr_1,c_qr_2) -> C_gg (ground diagonal).

    ground_id selects which diagonal element of the 4x4 matrix to return.
    """
    base = Path(__file__).parent
    pid = port_id
    suffix = "_sensing" if alternative_geometry else ""
    file_name = f"sweeps/clockmon_capacitance_library_{pid}{suffix}{_DEFAULT_SUFFIX}"
    file_path = base / file_name

    _, CMatrix, _ = _load_4x4_csv(str(file_path))
    wg_lib = waveguide_library()
    CMatrix = CMatrix.copy()
    CMatrix[:, 0, 0] -= wg_lib(deembed)
    CMatrix[:, 1, 1] -= wg_lib(deembed)

    c_qr_1, c_qr_2 = get_cqr(CMatrix)
    library = interpolate.LinearNDInterpolator((c_qr_1, c_qr_2), CMatrix[:, ground_id, ground_id])
    return library
