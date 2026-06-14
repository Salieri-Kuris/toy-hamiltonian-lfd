"""Plot the zero-field spectrum of the D3 two-band Hamiltonian."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from twse2_continuum_dos.d3_two_band import (
    D3TwoBandParameters,
    candidate_report,
    eigvals,
    line_spectrum,
    slope_reversal_candidate_parameters,
    strong_bending_candidate_parameters,
)


def band_mesh(params: D3TwoBandParameters, nk: int = 241) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return k-grid and lower/upper zero-field band energies."""
    k = np.linspace(-np.pi, np.pi, nk)
    lower = np.empty((nk, nk), dtype=float)
    upper = np.empty((nk, nk), dtype=float)
    for iy, k2 in enumerate(k):
        for ix, k1 in enumerate(k):
            lower[iy, ix], upper[iy, ix] = eigvals(float(k1), float(k2), params)
    return k, k, lower, upper


def path_from_nodes(
    nodes: list[tuple[str, np.ndarray]],
    points_per_segment: int = 420,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[tuple[str, float]]]:
    """Return a piecewise-linear path in (k1,k2)."""
    k_points: list[np.ndarray] = []
    x_values: list[float] = []
    ticks: list[tuple[str, float]] = [(nodes[0][0], 0.0)]
    x_current = 0.0
    for index, ((_, start), (label, stop)) in enumerate(zip(nodes[:-1], nodes[1:])):
        segment = np.linspace(0.0, 1.0, points_per_segment, endpoint=index == len(nodes) - 2)
        delta = stop - start
        length = float(np.linalg.norm(delta))
        if index > 0:
            segment = segment[1:]
        for value in segment:
            k_points.append(start + value * delta)
            x_values.append(x_current + value * length)
        x_current += length
        ticks.append((label, x_current))
    path = np.vstack(k_points)
    return path[:, 0], path[:, 1], np.asarray(x_values, dtype=float), ticks


def c2t_red_paths(points_per_segment: int = 420) -> list[tuple[str, np.ndarray, np.ndarray, np.ndarray, list[tuple[str, float]]]]:
    """Return the two red C2'T-invariant paths used for Dirac diagnostics."""
    gamma_path = path_from_nodes(
        [
            ("K'", np.array([-2.0 * np.pi / 3.0, 2.0 * np.pi / 3.0])),
            ("Gamma", np.array([0.0, 0.0])),
            ("K", np.array([2.0 * np.pi / 3.0, -2.0 * np.pi / 3.0])),
        ],
        points_per_segment=points_per_segment,
    )
    boundary_path = path_from_nodes(
        [
            ("K'", np.array([np.pi / 3.0, 2.0 * np.pi / 3.0])),
            ("M", np.array([np.pi / 2.0, np.pi / 2.0])),
            ("K", np.array([2.0 * np.pi / 3.0, np.pi / 3.0])),
        ],
        points_per_segment=points_per_segment,
    )
    return [
        ("C2'T red line: K'-Gamma-K", *gamma_path),
        ("C2'T boundary: K'-M-K", *boundary_path),
    ]


def selected_parameters(name: str) -> D3TwoBandParameters:
    if name == "strong":
        return strong_bending_candidate_parameters()
    if name == "slope_reversal":
        return slope_reversal_candidate_parameters()
    raise ValueError(f"unknown candidate: {name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", choices=["strong", "slope_reversal"], default="slope_reversal")
    args = parser.parse_args()

    params = selected_parameters(args.candidate)
    report = candidate_report(params, nk=151)
    selected = report["selected_dirac"]
    saddles = report["upper_saddles_above_dirac"]
    nearest_saddle = saddles[0] if saddles else None

    red_paths = c2t_red_paths(points_per_segment=480)
    k1_grid, k2_grid, lower, upper = band_mesh(params, nk=241)

    out_dir = ROOT / "figures"
    out_dir.mkdir(exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_path = out_dir / f"d3_two_band_zero_field_spectrum_{stamp}.png"

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.3), constrained_layout=True)

    for ax, (title, path_k1, path_k2, path_x, path_ticks) in zip(axes[0], red_paths):
        path_bands = line_spectrum(path_k1, path_k2, params)
        ax.plot(path_x, path_bands[:, 0], color="tab:blue", lw=1.1, label="lower")
        ax.plot(path_x, path_bands[:, 1], color="tab:orange", lw=1.1, label="upper")
        ax.axhline(selected["energy"], color="black", lw=0.8, alpha=0.6)
        if nearest_saddle is not None:
            ax.axhline(nearest_saddle["energy"], color="purple", lw=0.8, ls=":", alpha=0.8)
        for label, xpos in path_ticks:
            ax.axvline(xpos, color="black", lw=0.6, alpha=0.45)
        ax.set_xticks([xpos for _, xpos in path_ticks])
        ax.set_xticklabels([label for label, _ in path_ticks])
        ax.set_title(title)
        ax.set_xlabel("path")
        ax.set_ylabel("energy")
        ax.legend(fontsize=8, loc="best")

    common_levels = np.linspace(float(lower.min()), float(upper.max()), 32)
    im0 = axes[1, 0].imshow(
        lower,
        origin="lower",
        extent=[k1_grid[0], k1_grid[-1], k2_grid[0], k2_grid[-1]],
        cmap="viridis",
        aspect="equal",
    )
    axes[1, 0].contour(k1_grid, k2_grid, lower, levels=[selected["energy"]], colors="white", linewidths=0.7)
    axes[1, 0].plot([-np.pi, np.pi], [np.pi, -np.pi], color="red", lw=0.9, alpha=0.75)
    axes[1, 0].plot([0.0, np.pi], [np.pi, 0.0], color="red", lw=0.9, alpha=0.75)
    axes[1, 0].set_title("Lower band E-")
    axes[1, 0].set_xlabel("k1")
    axes[1, 0].set_ylabel("k2")
    fig.colorbar(im0, ax=axes[1, 0], label="energy")

    im1 = axes[1, 1].imshow(
        upper,
        origin="lower",
        extent=[k1_grid[0], k1_grid[-1], k2_grid[0], k2_grid[-1]],
        cmap="magma",
        aspect="equal",
    )
    axes[1, 1].contour(k1_grid, k2_grid, upper, levels=common_levels, colors="white", linewidths=0.25, alpha=0.32)
    axes[1, 1].contour(k1_grid, k2_grid, upper, levels=[selected["energy"]], colors="cyan", linewidths=0.8)
    axes[1, 1].plot([-np.pi, np.pi], [np.pi, -np.pi], color="red", lw=0.9, alpha=0.75)
    axes[1, 1].plot([0.0, np.pi], [np.pi, 0.0], color="red", lw=0.9, alpha=0.75)
    if nearest_saddle is not None:
        axes[1, 1].scatter([nearest_saddle["k1"]], [nearest_saddle["k2"]], s=42, color="cyan", edgecolor="black", lw=0.4)
        axes[1, 1].contour(k1_grid, k2_grid, upper, levels=[nearest_saddle["energy"]], colors="cyan", linewidths=0.8)
    axes[1, 1].set_title("Upper band E+")
    axes[1, 1].set_xlabel("k1")
    axes[1, 1].set_ylabel("k2")
    fig.colorbar(im1, ax=axes[1, 1], label="energy")

    fig.suptitle(
        f"D3 two-band zero-field spectrum ({args.candidate}), E_D={selected['energy']:.5f}, "
        f"nearest upper saddle dE={nearest_saddle['energy_minus_dirac']:.5f}" if nearest_saddle else "D3 two-band zero-field spectrum",
        fontsize=12,
    )
    fig.savefig(out_path, dpi=220)
    plt.close(fig)

    print(out_path.resolve())
    print(f"candidate={args.candidate}")
    print(f"E_D={selected['energy']:.9f}")
    print(f"Dirac roots={[(round(item['k1'], 9), round(item['k2'], 9)) for item in report['line_dirac_roots']]}")
    if nearest_saddle is not None:
        print(f"nearest_upper_saddle_dE={nearest_saddle['energy_minus_dirac']:.9f}")
        print(f"nearest_upper_saddle_k=({nearest_saddle['k1']:.9f}, {nearest_saddle['k2']:.9f})")


if __name__ == "__main__":
    main()
