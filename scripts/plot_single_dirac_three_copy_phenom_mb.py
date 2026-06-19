"""Plot the single-Dirac three-copy phenomenological MB model.

This follows ``single_dirac_three_copy_phenom_mb_manual.md``: the Landau
levels are semiclassical single-cone levels, while magnetic breakdown is a
phenomenological three-copy ring coupling in orbit space. It is not an exact
microscopic minimal-coupling Hamiltonian.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "hamiltonian_lfd_package_useful_docs_20260614-182601"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from twse2_continuum_dos.dos import FanResult, histogram_gaussian_dos, regrid_dos_to_density  # noqa: E402
from twse2_continuum_dos.holomorphic_cubic import row_normalized_contrast  # noqa: E402
from twse2_continuum_dos.phenomenology import (  # noqa: E402
    RouteABCResult,
    RouteABParams,
    RouteCParams,
    apply_route_abc_phenomenology,
)


@dataclass(frozen=True)
class PhenomParams:
    ED: float = 0.0
    v: float = 1.0
    alpha: float = -0.25
    K0: float = 1.0
    Nmax: int = 40
    Gamma_LL: float = 0.03
    t_scale: float = 0.06
    phi_MB: float = 0.5
    C_particle: float = 1.0
    C_hole: float = 12.0
    Ebeta: float = 0.45
    Gbeta: float = 0.10
    Abeta: float = 1.5

    def as_dict(self) -> dict[str, float | int]:
        out = asdict(self)
        out["Nmax"] = int(self.Nmax)
        return out


@dataclass(frozen=True)
class DOSData:
    energy_grid: np.ndarray
    b_values: np.ndarray
    n_grid: np.ndarray
    source_dos_EB: np.ndarray
    zero_field_dos: np.ndarray
    fan: FanResult
    route_abc: RouteABCResult


def positive_roots_quadratic(a: float, b: float, c: float, tol: float = 1e-12) -> list[float]:
    if abs(a) < tol:
        if abs(b) < tol:
            return []
        root = -c / b
        return [float(root)] if root >= -tol else []
    disc = b * b - 4.0 * a * c
    if disc < -tol:
        return []
    disc = max(float(disc), 0.0)
    roots = [(-b + np.sqrt(disc)) / (2.0 * a), (-b - np.sqrt(disc)) / (2.0 * a)]
    return sorted(float(root) for root in roots if root >= -tol)


def side_sign(side: str) -> int:
    if side == "particle":
        return 1
    if side == "hole":
        return -1
    raise ValueError("side must be 'particle' or 'hole'")


def orbit_radii(E: float, side: str, params: PhenomParams) -> list[float]:
    s = side_sign(side)
    return positive_roots_quadratic(params.alpha, s * params.v, -(float(E) - params.ED))


def choose_alpha_radius(E: float, side: str, params: PhenomParams) -> float:
    roots = orbit_radii(E, side, params)
    return float(roots[0]) if roots else float("nan")


def delta_k(E: float, side: str, params: PhenomParams) -> float:
    radius = choose_alpha_radius(E, side, params)
    if not np.isfinite(radius):
        return float("nan")
    return max(float(np.sqrt(3.0) * params.K0 - 2.0 * radius), 0.0)


def mb_probability(E: float, B: float, side: str, params: PhenomParams) -> float:
    if B <= 0.0:
        return 0.0
    Cgeom = params.C_particle if side == "particle" else params.C_hole
    dk = delta_k(E, side, params)
    if not np.isfinite(dk):
        return 0.0
    return float(np.exp(-float(Cgeom) * dk * dk / float(B)))


def dirac_LL_energy(N: int, B: float, side: str, params: PhenomParams) -> float:
    rN = np.sqrt(2.0 * float(B) * int(N))
    return float(params.ED + params.alpha * rN * rN + side_sign(side) * params.v * rN)


def split_three_copy_levels(B: float, side: str, params: PhenomParams) -> tuple[np.ndarray, np.ndarray]:
    levels: list[float] = [float(params.ED)] * 3
    weights: list[float] = [1.0] * 3
    for N in range(1, int(params.Nmax) + 1):
        E0 = dirac_LL_energy(N, B, side, params)
        probability = mb_probability(E0, B, side, params)
        amplitude = params.t_scale * np.sqrt(probability)
        for m in range(3):
            phase = (params.phi_MB + 2.0 * np.pi * m) / 3.0
            levels.append(float(E0 + 2.0 * amplitude * np.cos(phase)))
            weights.append(1.0)
    return np.asarray(levels, dtype=float), np.asarray(weights, dtype=float)


def beta_background_dos(energy_grid: np.ndarray, B: float, params: PhenomParams) -> np.ndarray:
    grid = np.asarray(energy_grid, dtype=float)
    lorentzian = (1.0 / np.pi) * params.Gbeta / ((grid - params.Ebeta) ** 2 + params.Gbeta**2)
    gate = 0.5 * (1.0 + np.tanh((grid - params.ED) / max(params.Gbeta, 1e-9)))
    probability = mb_probability(params.Ebeta, B, "particle", params)
    return params.Abeta * probability * lorentzian * gate


def zero_field_dos_profile(energy_grid: np.ndarray, params: PhenomParams, q_max: float = 3.0, n_q: int = 24000) -> np.ndarray:
    q = np.linspace(0.0, float(q_max), int(n_q))
    energies = []
    weights = []
    radial_weight = q.copy()
    for side in ("particle", "hole"):
        s = side_sign(side)
        branch = params.ED + params.alpha * q * q + s * params.v * q
        energies.append(branch)
        weights.append(3.0 * radial_weight)
    values = np.concatenate(energies)
    state_weights = np.concatenate(weights)
    profile = histogram_gaussian_dos(
        energy_grid,
        values,
        eta=max(float(params.Gamma_LL), 0.5 * float(params.Gbeta)),
        weights=state_weights,
        degeneracy=1.0,
    )
    profile += beta_background_dos(energy_grid, B=0.75, params=params)
    return np.maximum(profile, 0.0)


def source_dos_row(energy_grid: np.ndarray, B: float, params: PhenomParams) -> np.ndarray:
    dos = np.zeros_like(np.asarray(energy_grid, dtype=float))
    for side in ("particle", "hole"):
        levels, weights = split_three_copy_levels(B, side, params)
        degeneracy = float(B) / (2.0 * np.pi)
        dos += histogram_gaussian_dos(
            energy_grid,
            levels,
            eta=float(params.Gamma_LL),
            weights=weights,
            degeneracy=degeneracy,
        )
    dos += beta_background_dos(energy_grid, B, params)
    return np.maximum(dos, 0.0)


def build_dos_data(
    params: PhenomParams,
    *,
    B_values: np.ndarray,
    energy_grid: np.ndarray,
    n_grid: np.ndarray,
) -> DOSData:
    fields = np.asarray(B_values, dtype=float)
    energy = np.asarray(energy_grid, dtype=float)
    density = np.asarray(n_grid, dtype=float)
    source = np.vstack([source_dos_row(energy, float(B), params) for B in fields])
    zero_field = zero_field_dos_profile(energy, params)
    fan = regrid_dos_to_density(energy, fields, source, density, reference_energy=params.ED, out_of_range="nan")
    route_abc = apply_route_abc_phenomenology(
        energy_grid=energy,
        b_values=fields,
        source_dos_EB=source,
        zero_field_dos=zero_field,
        n_grid=density,
        reference_energy=params.ED,
        route_ab=RouteABParams(gamma0=0.006, a_sqrtB=0.004, lambda_D=0.012, zero_field_smooth_energy=0.015),
        route_c=RouteCParams(smooth_E=0.06, smooth_B=0.08, a_min=0.22, alpha=1.25),
    )
    return DOSData(
        energy_grid=energy,
        b_values=fields,
        n_grid=density,
        source_dos_EB=source,
        zero_field_dos=zero_field,
        fan=fan,
        route_abc=route_abc,
    )


def mb_diagnostics(params: PhenomParams, *, energy: float, B: float) -> dict[str, float]:
    p_particle = mb_probability(float(energy), float(B), "particle", params)
    p_hole = mb_probability(-abs(float(energy)), float(B), "hole", params)
    return {
        "energy_particle": float(energy),
        "energy_hole": -abs(float(energy)),
        "B": float(B),
        "P_particle": float(p_particle),
        "P_hole": float(p_hole),
        "particle_split": float(2.0 * params.t_scale * np.sqrt(p_particle)),
        "hole_split": float(2.0 * params.t_scale * np.sqrt(p_hole)),
        "delta_k_particle": float(delta_k(float(energy), "particle", params)),
        "delta_k_hole": float(delta_k(-abs(float(energy)), "hole", params)),
    }


def spectrum_points(B_values: np.ndarray, params: PhenomParams) -> tuple[np.ndarray, np.ndarray]:
    fields = []
    energies = []
    for B in np.asarray(B_values, dtype=float):
        for side in ("particle", "hole"):
            levels, _weights = split_three_copy_levels(float(B), side, params)
            fields.extend([float(B)] * len(levels))
            energies.extend(float(level) for level in levels)
    return np.asarray(fields, dtype=float), np.asarray(energies, dtype=float)


def robust_limits(values: np.ndarray, low: float = 1.0, high: float = 99.3) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 0.0, 1.0
    return float(np.percentile(finite, low)), float(np.percentile(finite, high))


def plot_spectrum(B_values: np.ndarray, params: PhenomParams, path: Path, energy_window: float) -> None:
    fields, energies = spectrum_points(B_values, params)
    fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    ax.scatter(fields, energies, s=4, c="black", alpha=0.58, linewidths=0)
    ax.axhline(params.ED, color="tab:red", lw=0.9, alpha=0.45)
    ax.axhline(params.Ebeta, color="tab:blue", lw=0.9, ls="--", alpha=0.50)
    ax.set_xlim(float(np.min(B_values)), float(np.max(B_values)))
    ax.set_ylim(-float(energy_window), float(energy_window))
    ax.set_xlabel("B")
    ax.set_ylabel("Energy E")
    ax.set_title("Single Dirac three-copy phenomenological MB spectrum")
    ax.grid(True, lw=0.4, alpha=0.25)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def plot_dos(data: DOSData, side_by_side_path: Path, dos_EB_path: Path, dos_nB_path: Path) -> None:
    density_EB = data.route_abc.d_density_EB
    visible_nB = data.route_abc.d_vis_nB
    eb_vmin, eb_vmax = robust_limits(density_EB)
    nb_contrast = row_normalized_contrast(np.nan_to_num(visible_nB, nan=0.0))

    fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    mesh = ax.pcolormesh(data.energy_grid, data.b_values, density_EB, shading="auto", cmap="magma", vmin=eb_vmin, vmax=eb_vmax)
    ax.set_xlabel("Energy E")
    ax.set_ylabel("B")
    ax.set_title("Phenomenological MB DOS: D(E,B)")
    fig.colorbar(mesh, ax=ax, label="density DOS")
    fig.savefig(dos_EB_path, dpi=240)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    mesh = ax.pcolormesh(data.n_grid, data.b_values, nb_contrast, shading="auto", cmap="magma", vmin=0.0, vmax=1.0)
    ax.set_xlabel("Density n")
    ax.set_ylabel("B")
    ax.set_title("Phenomenological MB DOS contrast: D(n,B)")
    fig.colorbar(mesh, ax=ax, label="row-normalized visible DOS")
    fig.savefig(dos_nB_path, dpi=240)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    mesh0 = axes[0].pcolormesh(data.energy_grid, data.b_values, density_EB, shading="auto", cmap="magma", vmin=eb_vmin, vmax=eb_vmax)
    axes[0].set_xlabel("Energy E")
    axes[0].set_ylabel("B")
    axes[0].set_title("D(E,B)")
    fig.colorbar(mesh0, ax=axes[0], label="density DOS")
    mesh1 = axes[1].pcolormesh(data.n_grid, data.b_values, nb_contrast, shading="auto", cmap="magma", vmin=0.0, vmax=1.0)
    axes[1].set_xlabel("Density n")
    axes[1].set_ylabel("B")
    axes[1].set_title("D(n,B) visible contrast")
    fig.colorbar(mesh1, ax=axes[1], label="row-normalized visible DOS")
    fig.savefig(side_by_side_path, dpi=240)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--b-min", type=float, default=0.08)
    parser.add_argument("--b-max", type=float, default=1.20)
    parser.add_argument("--n-b", type=int, default=120)
    parser.add_argument("--energy-window", type=float, default=1.35)
    parser.add_argument("--n-energy", type=int, default=1201)
    parser.add_argument("--n-min", type=float, default=-0.28)
    parser.add_argument("--n-max", type=float, default=0.28)
    parser.add_argument("--n-density", type=int, default=801)
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "outputs" / "single_dirac_three_copy_phenom_mb")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    params = PhenomParams()
    B_values = np.linspace(float(args.b_min), float(args.b_max), int(args.n_b))
    energy_grid = np.linspace(-float(args.energy_window), float(args.energy_window), int(args.n_energy))
    n_grid = np.linspace(float(args.n_min), float(args.n_max), int(args.n_density))

    data = build_dos_data(params, B_values=B_values, energy_grid=energy_grid, n_grid=n_grid)

    stem = f"single_dirac_three_copy_phenom_mb_{timestamp}"
    npz_path = args.output_dir / f"{stem}.npz"
    summary_path = args.output_dir / f"{stem}.json"
    spectrum_path = args.output_dir / f"{stem}_spectrum.png"
    dos_EB_path = args.output_dir / f"{stem}_dos_EB.png"
    dos_nB_path = args.output_dir / f"{stem}_dos_nB.png"
    side_by_side_path = args.output_dir / f"{stem}_dos_EB_nB_side_by_side.png"

    np.savez(
        npz_path,
        B_values=data.b_values,
        energy_grid=data.energy_grid,
        n_grid=data.n_grid,
        source_dos_EB=data.source_dos_EB,
        zero_field_dos=data.zero_field_dos,
        d_density_EB=data.route_abc.d_density_EB,
        d_vis_EB=data.route_abc.d_vis_EB,
        d_density_nB=data.route_abc.d_density_nB,
        d_vis_nB=data.route_abc.d_vis_nB,
        mu_nB=data.route_abc.mu_nB,
        **params.as_dict(),
    )
    plot_spectrum(B_values, params, spectrum_path, float(args.energy_window))
    plot_dos(data, side_by_side_path, dos_EB_path, dos_nB_path)

    diagnostics = mb_diagnostics(params, energy=params.Ebeta, B=0.75)
    summary = {
        "created_at": timestamp,
        "manual": "single_dirac_three_copy_phenom_mb_manual.md",
        "model": "single Dirac cone, three C3 copies, phenomenological orbit-space MB",
        "note": "MB coupling is phenomenological in LL/copy space, not exact microscopic H0(Pi).",
        "params": params.as_dict(),
        "diagnostics_at_Ebeta_B0p75": diagnostics,
        "B_min": float(args.b_min),
        "B_max": float(args.b_max),
        "n_B": int(args.n_b),
        "energy_window": float(args.energy_window),
        "n_energy": int(args.n_energy),
        "n_min": float(args.n_min),
        "n_max": float(args.n_max),
        "n_density": int(args.n_density),
        "npz": str(npz_path),
        "summary": str(summary_path),
        "spectrum_figure": str(spectrum_path),
        "dos_EB_figure": str(dos_EB_path),
        "dos_nB_figure": str(dos_nB_path),
        "dos_side_by_side_figure": str(side_by_side_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
