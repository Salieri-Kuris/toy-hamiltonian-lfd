"""Physically constrained two-band D3 reservoir model.

The construction keeps one two-band block per valley,

    H(k) = d0(k) sigma_0 + M(k) sigma_z + G(k) sigma_x,

with all three terms represented by finite Fourier harmonics.  The default
places the Dirac point on a generic C2T-invariant line rather than at K/K',
which gives enough particle-hole asymmetry to keep the hole side clean while
leaving a single Gamma-centered particle-side reservoir.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.ndimage as ndi

from .d3_two_band import S0, SX, SZ
from .peierls_hofstadter import TBModel


@dataclass(frozen=True)
class PhysicalReservoirParams:
    """Low-dimensional knobs for a single particle-side reservoir."""

    root: float = 1.00
    reservoir_energy: float = 0.086
    reservoir_depth: float = 0.070
    alpha_velocity: float = 0.220
    offdiag_velocity: float = 0.500
    boundary_dispersion: float = 0.040
    side_raise: float = 0.0
    profile_power: float = 0.45


@dataclass(frozen=True)
class ReservoirWindowMetrics:
    """Zero-field phase-space diagnostics in the particle reservoir window."""

    center_window_area: float
    side_window_area: float
    max_side_shell_area: float
    center_to_side_ratio: float


@dataclass(frozen=True)
class ParticleBreakdownMetrics:
    """Zero-field diagnostics for particle-only magnetic-breakdown geometry."""

    hole_background_area: float
    hole_max_shell_area: float
    low_particle_alpha_components: int
    low_particle_background_components: int
    reservoir_components: int
    neck_energy: float
    particle_hole_asymmetry: float


def c1_harmonic(k1: np.ndarray | float, k2: np.ndarray | float) -> np.ndarray:
    """Return the first triangular D3 cosine harmonic."""
    k1a = np.asarray(k1)
    k2a = np.asarray(k2)
    return np.cos(k1a) + np.cos(k2a) + np.cos(k1a + k2a)


def c2_harmonic(k1: np.ndarray | float, k2: np.ndarray | float) -> np.ndarray:
    """Return the second triangular D3 cosine harmonic."""
    k1a = np.asarray(k1)
    k2a = np.asarray(k2)
    return np.cos(2.0 * k1a) + np.cos(2.0 * k2a) + np.cos(2.0 * (k1a + k2a))


def g_harmonic(k1: np.ndarray | float, k2: np.ndarray | float) -> np.ndarray:
    """Return a C2T-line-vanishing off-diagonal harmonic."""
    return np.sin(np.asarray(k1) + np.asarray(k2))


def _root_c1(params: PhysicalReservoirParams) -> float:
    return float(c1_harmonic(params.root, -params.root))


def _root_c2(params: PhysicalReservoirParams) -> float:
    return float(c2_harmonic(params.root, -params.root))


def _mass_normalizer(params: PhysicalReservoirParams) -> float:
    # Along k2=-k1, c1=1+2 cos(k1).  This normalizes the Dirac slope so the
    # alpha_velocity knob has a direct local meaning near the chosen root.
    slope = abs(-2.0 * np.sin(params.root))
    return max(slope, 1.0e-12)


def _c2_normalizer(params: PhysicalReservoirParams) -> float:
    return max(abs(3.0 - _root_c2(params)), 1.0e-12)


def mass_value(
    k1: np.ndarray | float,
    k2: np.ndarray | float,
    params: PhysicalReservoirParams,
) -> np.ndarray:
    """Return M(k), with roots at (root, -root) and (-root, root)."""
    return params.alpha_velocity * (c1_harmonic(k1, k2) - _root_c1(params)) / _mass_normalizer(params)


def offdiag_value(
    k1: np.ndarray | float,
    k2: np.ndarray | float,
    params: PhysicalReservoirParams,
) -> np.ndarray:
    """Return G(k), which vanishes on the representative C2T line."""
    return params.offdiag_velocity * g_harmonic(k1, k2)


def target_upper_value(
    k1: np.ndarray | float,
    k2: np.ndarray | float,
    params: PhysicalReservoirParams,
) -> np.ndarray:
    """Return the upper branch implied by the smooth harmonic Hamiltonian."""
    d0 = d0_value(k1, k2, params)
    gap = np.sqrt(mass_value(k1, k2, params) ** 2 + offdiag_value(k1, k2, params) ** 2)
    return d0 + gap


def scalar_tilt(params: PhysicalReservoirParams) -> float:
    """Return the d0 slope that places the Gamma reservoir at the target energy."""
    x_gamma = (3.0 - _root_c1(params)) / _mass_normalizer(params)
    y_gamma = (3.0 - _root_c2(params)) / _c2_normalizer(params)
    if x_gamma <= 0.0:
        raise ValueError("chosen root does not put Gamma on the particle side")
    tilt = (params.reservoir_energy - params.boundary_dispersion * y_gamma) / x_gamma - params.alpha_velocity
    limit = 0.98 * params.alpha_velocity
    return float(np.clip(tilt, -limit, limit))


def d0_value(
    k1: np.ndarray | float,
    k2: np.ndarray | float,
    params: PhysicalReservoirParams,
) -> np.ndarray:
    """Return a smooth scalar tilt that keeps the Dirac energy at half filling."""
    x = (c1_harmonic(k1, k2) - _root_c1(params)) / _mass_normalizer(params)
    y = (c2_harmonic(k1, k2) - _root_c2(params)) / _c2_normalizer(params)
    return scalar_tilt(params) * x + params.boundary_dispersion * y


def hamiltonian(
    k1: float,
    k2: float,
    params: PhysicalReservoirParams,
) -> np.ndarray:
    """Return the 2 x 2 Bloch Hamiltonian for one valley."""
    d0 = float(d0_value(k1, k2, params))
    mass = float(mass_value(k1, k2, params))
    offdiag = float(offdiag_value(k1, k2, params))
    return np.array([[d0 + mass, offdiag], [offdiag, d0 - mass]], dtype=float)


def bands(
    k1: np.ndarray | float,
    k2: np.ndarray | float,
    params: PhysicalReservoirParams,
) -> tuple[np.ndarray, np.ndarray]:
    """Return lower and upper eigenvalues."""
    d0 = d0_value(k1, k2, params)
    gap = np.sqrt(mass_value(k1, k2, params) ** 2 + offdiag_value(k1, k2, params) ** 2)
    return d0 - gap, d0 + gap


def construct_single_reservoir_candidate(root: float = 1.00) -> PhysicalReservoirParams:
    """Return a deterministic candidate with one dominant particle reservoir."""
    return PhysicalReservoirParams(root=float(root))


def reservoir_window_metrics(
    params: PhysicalReservoirParams,
    *,
    nk: int = 181,
    energy_window: tuple[float, float] = (0.036, 0.088),
    center_radius: float = 1.18,
    dirac_exclusion_radius: float = 0.55,
    shell_levels: tuple[float, ...] = (0.042, 0.050, 0.058, 0.070, 0.082),
    shell_width: float = 0.006,
) -> ReservoirWindowMetrics:
    """Measure whether the active particle window has only one central family."""
    grid = np.linspace(-np.pi, np.pi, int(nk), endpoint=False)
    k1, k2 = np.meshgrid(grid, grid, indexing="ij")
    _, upper = bands(k1, k2, params)
    in_window = (upper >= energy_window[0]) & (upper <= energy_window[1])
    dirac_distance = np.minimum(
        np.hypot(k1 - params.root, k2 + params.root),
        np.hypot(k1 + params.root, k2 - params.root),
    )
    alpha_neighborhood = dirac_distance < dirac_exclusion_radius
    non_alpha_window = in_window & ~alpha_neighborhood
    labels, count = ndi.label(non_alpha_window, structure=np.ones((3, 3), dtype=bool))
    gamma_index = (int(np.argmin(np.abs(grid))), int(np.argmin(np.abs(grid))))
    gamma_label = int(labels[gamma_index])
    if gamma_label == 0:
        center = np.zeros_like(non_alpha_window, dtype=bool)
    else:
        center = labels == gamma_label
    side = non_alpha_window & ~center
    side_shell_areas = [
        float(np.mean((np.abs(upper - level) <= shell_width) & side))
        for level in shell_levels
    ]
    center_area = float(np.mean(center))
    side_area = float(np.mean(side))
    return ReservoirWindowMetrics(
        center_window_area=center_area,
        side_window_area=side_area,
        max_side_shell_area=float(max(side_shell_areas, default=0.0)),
        center_to_side_ratio=center_area / max(side_area, 1.0e-12),
    )


def _component_count(mask: np.ndarray, *, min_pixels: int = 10) -> int:
    labels, count = ndi.label(mask, structure=np.ones((3, 3), dtype=bool))
    kept = 0
    for label in range(1, count + 1):
        if int(np.count_nonzero(labels == label)) >= int(min_pixels):
            kept += 1
    return kept


def particle_breakdown_metrics(
    params: PhysicalReservoirParams,
    *,
    nk: int = 201,
    low_particle_level: float = 0.014,
    reservoir_window: tuple[float, float] = (0.046, 0.088),
    shell_width: float = 0.004,
    dirac_radius: float = 0.42,
    hole_window: tuple[float, float] = (-0.085, -0.020),
    hole_shell_levels: tuple[float, ...] = (-0.080, -0.060, -0.040, -0.025),
) -> ParticleBreakdownMetrics:
    """Measure whether only the particle side has alpha/reservoir breakdown geometry."""
    grid = np.linspace(-np.pi, np.pi, int(nk), endpoint=False)
    k1, k2 = np.meshgrid(grid, grid, indexing="ij")
    lower, upper = bands(k1, k2, params)
    dirac_distance = np.minimum(
        np.hypot(k1 - params.root, k2 + params.root),
        np.hypot(k1 + params.root, k2 - params.root),
    )
    alpha_neighborhood = dirac_distance < float(dirac_radius)

    particle_shell = np.abs(upper - float(low_particle_level)) <= float(shell_width)
    low_alpha_components = _component_count(particle_shell & alpha_neighborhood, min_pixels=6)
    low_background_components = _component_count(particle_shell & ~alpha_neighborhood, min_pixels=10)

    reservoir_mask = (
        (upper >= reservoir_window[0])
        & (upper <= reservoir_window[1])
        & (np.hypot(k1, k2) < 1.35)
        & ~alpha_neighborhood
    )
    reservoir_components = _component_count(reservoir_mask, min_pixels=20)

    shell_levels = np.linspace(low_particle_level, reservoir_window[1], 24)
    background_areas = [
        float(np.mean((np.abs(upper - level) <= shell_width) & ~alpha_neighborhood))
        for level in shell_levels
    ]
    onset_indices = [index for index, area in enumerate(background_areas) if area > 1.0e-4]
    neck_energy = float(shell_levels[onset_indices[0]]) if onset_indices else float("nan")

    hole_mask = (lower >= hole_window[0]) & (lower <= hole_window[1]) & ~alpha_neighborhood
    hole_shell_areas = [
        float(np.mean((np.abs(lower - level) <= shell_width) & ~alpha_neighborhood))
        for level in hole_shell_levels
    ]

    particle_area = float(np.mean((upper >= low_particle_level) & (upper <= reservoir_window[1])))
    hole_area = float(np.mean(hole_mask))
    return ParticleBreakdownMetrics(
        hole_background_area=hole_area,
        hole_max_shell_area=float(max(hole_shell_areas, default=0.0)),
        low_particle_alpha_components=low_alpha_components,
        low_particle_background_components=low_background_components,
        reservoir_components=reservoir_components,
        neck_energy=neck_energy,
        particle_hole_asymmetry=particle_area / max(hole_area, 1.0e-12),
    )


def tb_model(params: PhysicalReservoirParams, *, atol: float = 1.0e-14) -> TBModel:
    """Return the exact finite-hopping model for the smooth harmonic ansatz."""
    c1_coeffs = {
        (1, 0): 0.5,
        (-1, 0): 0.5,
        (0, 1): 0.5,
        (0, -1): 0.5,
        (1, 1): 0.5,
        (-1, -1): 0.5,
    }
    sin_s_coeffs = {(1, 1): 1.0 / (2.0j), (-1, -1): -1.0 / (2.0j)}
    c2_coeffs = {
        (2, 0): 0.5,
        (-2, 0): 0.5,
        (0, 2): 0.5,
        (0, -2): 0.5,
        (2, 2): 0.5,
        (-2, -2): 0.5,
    }
    const_key = (0, 0)
    norm = _mass_normalizer(params)
    c2_norm = _c2_normalizer(params)
    root_c1 = _root_c1(params)
    root_c2 = _root_c2(params)
    tilt = scalar_tilt(params)

    coeffs: dict[tuple[int, int], np.ndarray] = {}

    def add(key: tuple[int, int], mat: np.ndarray) -> None:
        coeffs[key] = coeffs.get(key, np.zeros((2, 2), dtype=complex)) + mat

    for key, value in c1_coeffs.items():
        add(key, (tilt / norm) * value * S0)
        add(key, (params.alpha_velocity / norm) * value * SZ)
    add(const_key, -(tilt * root_c1 / norm) * S0)
    add(const_key, -(params.alpha_velocity * root_c1 / norm) * SZ)
    for key, value in c2_coeffs.items():
        add(key, (params.boundary_dispersion / c2_norm) * value * S0)
    add(const_key, -(params.boundary_dispersion * root_c2 / c2_norm) * S0)
    for key, value in sin_s_coeffs.items():
        add(key, params.offdiag_velocity * value * SX)

    hoppings = [
        (dm, dn, np.asarray(mat, dtype=complex))
        for (dm, dn), mat in sorted(coeffs.items())
        if np.linalg.norm(mat) > atol
    ]
    return TBModel(n_orb=2, hoppings=hoppings)
