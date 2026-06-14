"""Zero-field high-symmetry and contour diagnostics for D3 two-band models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.ndimage as ndi

from .d3_two_band import (
    D3TwoBandParameters,
    candidate_report,
    d0_term,
    g_harmonic,
    line_spectrum,
    m_term,
)


@dataclass(frozen=True)
class HighSymmetryPathSpectrum:
    """Band energies along one named high-symmetry path."""

    name: str
    k1: np.ndarray
    k2: np.ndarray
    x: np.ndarray
    ticks: list[tuple[str, float]]
    energies: np.ndarray


@dataclass(frozen=True)
class HighSymmetrySpectrum:
    """High-symmetry spectra and reference zero-field energies."""

    paths: list[HighSymmetryPathSpectrum]
    energy_dirac: float
    upper_saddle_energy: float | None


@dataclass(frozen=True)
class BandMesh:
    """Regular zero-field upper-band mesh."""

    k: np.ndarray
    k1: np.ndarray
    k2: np.ndarray
    upper: np.ndarray


@dataclass(frozen=True)
class ContourComponent:
    """One connected component of a narrow equal-energy shell."""

    label: int
    pixel_count: int
    centroid_k1: float
    centroid_k2: float
    boundary_touch: bool
    c3_self_overlap: float
    beta_like: bool
    alpha_like: bool


@dataclass(frozen=True)
class ContourLevelDiagnostic:
    """Contour component diagnostic at one particle-side energy level."""

    level: float
    offset: float
    component_count: int
    alpha_candidate_count: int
    beta_candidate_count: int
    alpha_mask_fraction: float
    beta_mask_fraction: float
    components: list[ContourComponent]


@dataclass(frozen=True)
class ContourDiagnosticResult:
    """Particle-side contour and alpha/beta mask diagnostic result."""

    energy_dirac: float
    levels: np.ndarray
    by_level: list[ContourLevelDiagnostic]
    well_defined_alpha: bool


def _path_from_nodes(
    nodes: list[tuple[str, np.ndarray]],
    points_per_segment: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[tuple[str, float]]]:
    points: list[np.ndarray] = []
    x_values: list[float] = []
    ticks: list[tuple[str, float]] = [(nodes[0][0], 0.0)]
    x_current = 0.0
    for index, ((_, start), (label, stop)) in enumerate(zip(nodes[:-1], nodes[1:])):
        segment = np.linspace(0.0, 1.0, points_per_segment, endpoint=index == len(nodes) - 2)
        if index > 0:
            segment = segment[1:]
        delta = stop - start
        length = float(np.linalg.norm(delta))
        for value in segment:
            points.append(start + value * delta)
            x_values.append(x_current + value * length)
        x_current += length
        ticks.append((label, x_current))
    path = np.vstack(points)
    return path[:, 0], path[:, 1], np.asarray(x_values, dtype=float), ticks


def high_symmetry_spectrum(
    params: D3TwoBandParameters,
    *,
    points_per_segment: int = 360,
) -> HighSymmetrySpectrum:
    """Return spectra on the two red C2T lines used by the D3 workflow."""
    root = 2.0 * np.pi / 3.0
    nodes = [
        (
            "K-Gamma-K'",
            [
                ("K", np.array([root, -root])),
                ("Gamma", np.array([0.0, 0.0])),
                ("K'", np.array([-root, root])),
            ],
        ),
        (
            "K'-M-K",
            [
                ("K'", np.array([np.pi / 3.0, 2.0 * np.pi / 3.0])),
                ("M", np.array([np.pi / 2.0, np.pi / 2.0])),
                ("K", np.array([2.0 * np.pi / 3.0, np.pi / 3.0])),
            ],
        ),
    ]
    paths: list[HighSymmetryPathSpectrum] = []
    for name, path_nodes in nodes:
        k1, k2, x, ticks = _path_from_nodes(path_nodes, points_per_segment)
        paths.append(
            HighSymmetryPathSpectrum(
                name=name,
                k1=k1,
                k2=k2,
                x=x,
                ticks=ticks,
                energies=line_spectrum(k1, k2, params),
            )
        )

    report = candidate_report(params, nk=101)
    saddles = report["upper_saddles_above_dirac"]
    return HighSymmetrySpectrum(
        paths=paths,
        energy_dirac=float(report["selected_dirac"]["energy"]),
        upper_saddle_energy=float(saddles[0]["energy"]) if saddles else None,
    )


def upper_band_mesh(params: D3TwoBandParameters, *, nk: int = 301) -> BandMesh:
    """Return a vectorized upper-band mesh on the periodic square cell."""
    k = np.linspace(-np.pi, np.pi, int(nk), endpoint=False)
    k1, k2 = np.meshgrid(k, k, indexing="xy")
    gap = np.sqrt(m_term(k1, k2, params) ** 2 + (params.lam * g_harmonic(k1, k2)) ** 2)
    upper = d0_term(k1, k2, params) + gap
    return BandMesh(k=k, k1=k1, k2=k2, upper=upper)


def _angle_mean(values: np.ndarray) -> float:
    z = np.mean(np.exp(1.0j * values))
    return float(np.angle(z))


def _wrap_index(index: np.ndarray, size: int) -> np.ndarray:
    return np.mod(index, size).astype(int)


def _c3_index_overlap(mask: np.ndarray, component: np.ndarray) -> float:
    n = mask.shape[0]
    rows, cols = np.nonzero(component)
    if rows.size == 0:
        return 0.0
    step = 2.0 * np.pi / n
    k1 = -np.pi + step * cols
    k2 = -np.pi + step * rows
    k1_rot = k2
    k2_rot = -k1 - k2
    col_rot = _wrap_index(np.rint((k1_rot + np.pi) / step).astype(int), n)
    row_rot = _wrap_index(np.rint((k2_rot + np.pi) / step).astype(int), n)
    hits = component[row_rot, col_rot]
    return float(np.count_nonzero(hits) / rows.size)


def _component_summary(
    labels: np.ndarray,
    label: int,
    mesh: BandMesh,
    *,
    total_pixels: int,
    c3_alpha_threshold: float,
    max_alpha_fraction: float,
) -> ContourComponent:
    component = labels == int(label)
    rows, cols = np.nonzero(component)
    pixel_count = int(rows.size)
    centroid_k1 = _angle_mean(mesh.k1[component])
    centroid_k2 = _angle_mean(mesh.k2[component])
    boundary_touch = bool(
        np.any(rows == 0)
        or np.any(cols == 0)
        or np.any(rows == labels.shape[0] - 1)
        or np.any(cols == labels.shape[1] - 1)
    )
    c3_overlap = _c3_index_overlap(labels > 0, component)
    fraction = pixel_count / float(total_pixels)
    alpha_like = bool((not boundary_touch) and c3_overlap >= c3_alpha_threshold and fraction <= max_alpha_fraction)
    beta_like = bool(boundary_touch or not alpha_like)
    return ContourComponent(
        label=int(label),
        pixel_count=pixel_count,
        centroid_k1=centroid_k1,
        centroid_k2=centroid_k2,
        boundary_touch=boundary_touch,
        c3_self_overlap=c3_overlap,
        beta_like=beta_like,
        alpha_like=alpha_like,
    )


def contour_component_diagnostics(
    params: D3TwoBandParameters,
    *,
    offsets: np.ndarray | list[float] | tuple[float, ...] = (0.018, 0.026, 0.034, 0.046),
    nk: int = 301,
    shell_width: float = 0.004,
    min_pixels: int = 12,
    c3_alpha_threshold: float = 0.55,
    max_alpha_fraction: float = 0.09,
) -> ContourDiagnosticResult:
    """Diagnose particle-side contour components and a coarse alpha/beta split.

    The mask is a narrow equal-energy shell diagnostic. It is not a strict
    Hilbert-space projector; it is meant to decide whether the zero-field phase
    space has a unique, stable alpha component before expensive Hofstadter runs.
    """
    report = candidate_report(params, nk=101)
    e_dirac = float(report["selected_dirac"]["energy"])
    levels = e_dirac + np.asarray(offsets, dtype=float)
    mesh = upper_band_mesh(params, nk=int(nk))
    by_level: list[ContourLevelDiagnostic] = []
    total_pixels = int(mesh.upper.size)
    for level, offset in zip(levels, np.asarray(offsets, dtype=float)):
        shell = np.abs(mesh.upper - float(level)) <= float(shell_width)
        labels, count = ndi.label(shell, structure=np.ones((3, 3), dtype=int))
        components: list[ContourComponent] = []
        for label in range(1, count + 1):
            if int(np.count_nonzero(labels == label)) < int(min_pixels):
                labels[labels == label] = 0
                continue
            components.append(
                _component_summary(
                    labels,
                    label,
                    mesh,
                    total_pixels=total_pixels,
                    c3_alpha_threshold=float(c3_alpha_threshold),
                    max_alpha_fraction=float(max_alpha_fraction),
                )
            )
        alpha_pixels = sum(component.pixel_count for component in components if component.alpha_like)
        beta_pixels = sum(component.pixel_count for component in components if component.beta_like)
        by_level.append(
            ContourLevelDiagnostic(
                level=float(level),
                offset=float(offset),
                component_count=len(components),
                alpha_candidate_count=sum(1 for component in components if component.alpha_like),
                beta_candidate_count=sum(1 for component in components if component.beta_like),
                alpha_mask_fraction=float(alpha_pixels / total_pixels),
                beta_mask_fraction=float(beta_pixels / total_pixels),
                components=components,
            )
        )

    well_defined = bool(by_level and all(item.alpha_candidate_count == 1 for item in by_level))
    return ContourDiagnosticResult(
        energy_dirac=e_dirac,
        levels=levels,
        by_level=by_level,
        well_defined_alpha=well_defined,
    )
