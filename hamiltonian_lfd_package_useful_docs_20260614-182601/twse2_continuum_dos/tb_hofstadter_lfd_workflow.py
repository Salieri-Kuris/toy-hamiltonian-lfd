"""Small workflow helpers for TB Peierls-Hofstadter LFD checks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .d3_two_band import (
    D3TwoBandParameters,
    slope_reversal_candidate_parameters,
    strong_bending_candidate_parameters,
    tb_model,
)
from .dos import FanResult, histogram_gaussian_dos, regrid_dos_to_density
from .peierls_hofstadter import magnetic_field_tesla_from_q, spectrum_at_flux


@dataclass(frozen=True)
class PureHofstadterDOSResult:
    """Pure Peierls-Hofstadter DOS and density-axis fan data."""

    q_values: np.ndarray
    b_values: np.ndarray
    energy_grid: np.ndarray
    dos_EB: np.ndarray
    fan: FanResult
    eta: float
    kmesh: tuple[int, int]


def selected_d3_parameters(name: str) -> D3TwoBandParameters:
    """Return a named D3 two-band candidate used by the TB manual workflow."""
    if name == "strong":
        return strong_bending_candidate_parameters()
    if name == "slope_reversal":
        return slope_reversal_candidate_parameters()
    raise ValueError(f"unknown D3 candidate: {name}")


def compute_pure_hofstadter_dos(
    params: D3TwoBandParameters,
    *,
    q_values: np.ndarray,
    energy_grid: np.ndarray,
    n_grid: np.ndarray,
    eta: float,
    kmesh: tuple[int, int],
    p: int = 1,
    reference_energy: float | None = None,
    degeneracy: float | None = None,
) -> PureHofstadterDOSResult:
    """Compute total DOS and total density fan from the same Peierls spectrum.

    This helper intentionally does not add a phenomenological beta reservoir,
    hand-drawn fan line, or magnetic-breakdown replica spectrum.
    """
    if eta <= 0.0:
        raise ValueError("eta must be positive")
    q_arr = np.asarray(q_values, dtype=int)
    if q_arr.ndim != 1 or q_arr.size == 0:
        raise ValueError("q_values must be a non-empty one-dimensional array")
    if np.any(q_arr <= 0):
        raise ValueError("q_values must be positive")

    fields = np.asarray([magnetic_field_tesla_from_q(int(q), p=p) for q in q_arr], dtype=float)
    order = np.argsort(fields)
    q_sorted = q_arr[order]
    b_sorted = fields[order]
    model = tb_model(params)
    total_degeneracy = float(params.internal_degeneracy if degeneracy is None else degeneracy)

    dos_rows = np.empty((q_sorted.size, len(energy_grid)), dtype=float)
    for index, q in enumerate(q_sorted):
        spectrum = spectrum_at_flux(
            model,
            p=int(p),
            q=int(q),
            kmesh=kmesh,
            method="dense_full",
            return_eigenvectors=False,
        )
        values = np.asarray(spectrum.eigenvalues, dtype=float).reshape(-1)
        weights = np.full(values.size, float(spectrum.state_weight_per_original_cell))
        dos_rows[index] = histogram_gaussian_dos(
            energy_grid,
            values,
            eta=float(eta),
            weights=weights,
            degeneracy=total_degeneracy,
        )

    fan = regrid_dos_to_density(energy_grid, b_sorted, dos_rows, n_grid, reference_energy=reference_energy)
    return PureHofstadterDOSResult(
        q_values=q_sorted,
        b_values=b_sorted,
        energy_grid=np.asarray(energy_grid, dtype=float),
        dos_EB=dos_rows,
        fan=fan,
        eta=float(eta),
        kmesh=tuple(int(item) for item in kmesh),
    )
