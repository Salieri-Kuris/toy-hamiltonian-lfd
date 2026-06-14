"""Plot zero-field and quick Hofstadter diagnostics for the physical reservoir model."""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from twse2_continuum_dos.d3_physical_reservoir import (
    bands,
    construct_single_reservoir_candidate,
    particle_breakdown_metrics,
    reservoir_window_metrics,
    tb_model,
)
from twse2_continuum_dos.dos import histogram_gaussian_dos, regrid_dos_to_density
from twse2_continuum_dos.peierls_hofstadter import magnetic_field_tesla_from_q, spectrum_at_flux


def path_from_nodes(
    nodes: list[tuple[str, np.ndarray]],
    points_per_segment: int = 480,
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
    return path[:, 0], path[:, 1], np.asarray(x_values), ticks


def high_symmetry_paths(root: float) -> list[tuple[str, np.ndarray, np.ndarray, np.ndarray, list[tuple[str, float]]]]:
    return [
        (
            "C2T line through shifted Dirac",
            *path_from_nodes(
                [
                    ("-D", np.array([-root, root])),
                    ("Gamma", np.array([0.0, 0.0])),
                    ("D", np.array([root, -root])),
                ]
            ),
        ),
        (
            "K'-M-K boundary",
            *path_from_nodes(
                [
                    ("K'", np.array([np.pi / 3.0, 2.0 * np.pi / 3.0])),
                    ("M", np.array([np.pi / 2.0, np.pi / 2.0])),
                    ("K", np.array([2.0 * np.pi / 3.0, np.pi / 3.0])),
                ]
            ),
        ),
    ]


def plot_zero_field(stamp: str) -> Path:
    params = construct_single_reservoir_candidate()
    figure_dir = ROOT / "figures"
    figure_dir.mkdir(exist_ok=True)
    out_path = figure_dir / f"d3_physical_reservoir_zero_field_{stamp}.png"

    grid = np.linspace(-np.pi, np.pi, 361, endpoint=False)
    k1, k2 = np.meshgrid(grid, grid, indexing="ij")
    lower, upper = bands(k1, k2, params)
    paths = high_symmetry_paths(params.root)

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.2), constrained_layout=True)
    for ax, (title, path_k1, path_k2, x, ticks) in zip(axes[0], paths):
        lo, up = bands(path_k1, path_k2, params)
        ax.plot(x, lo, color="tab:blue", lw=1.2, label="lower")
        ax.plot(x, up, color="tab:orange", lw=1.2, label="upper")
        ax.axhline(0.0, color="black", lw=0.8, alpha=0.65)
        ax.axhline(params.reservoir_energy, color="purple", lw=0.8, ls=":", alpha=0.9)
        for _, xpos in ticks:
            ax.axvline(xpos, color="black", lw=0.55, alpha=0.35)
        ax.set_xticks([xpos for _, xpos in ticks])
        ax.set_xticklabels([label for label, _ in ticks])
        ax.set_title(title)
        ax.set_ylabel("E - E_D")
        ax.legend(fontsize=8)

    im0 = axes[1, 0].imshow(
        lower.T,
        origin="lower",
        extent=[grid[0], grid[-1], grid[0], grid[-1]],
        cmap="viridis",
        aspect="equal",
    )
    axes[1, 0].contour(grid, grid, lower.T, levels=[-0.08, -0.06, -0.04, -0.02], colors="white", linewidths=0.35)
    axes[1, 0].scatter([params.root, -params.root], [-params.root, params.root], s=28, c="cyan", edgecolor="black", lw=0.4)
    axes[1, 0].plot([-np.pi, np.pi], [np.pi, -np.pi], color="red", lw=0.8, alpha=0.7)
    axes[1, 0].set_title("Lower band: clean hole side")
    axes[1, 0].set_xlabel("k1")
    axes[1, 0].set_ylabel("k2")
    fig.colorbar(im0, ax=axes[1, 0], label="E - E_D")

    im1 = axes[1, 1].imshow(
        upper.T,
        origin="lower",
        extent=[grid[0], grid[-1], grid[0], grid[-1]],
        cmap="magma",
        aspect="equal",
    )
    axes[1, 1].contour(grid, grid, upper.T, levels=[0.014, 0.043, 0.060, 0.080], colors="white", linewidths=0.55)
    axes[1, 1].scatter([params.root, -params.root], [-params.root, params.root], s=28, c="cyan", edgecolor="black", lw=0.4)
    axes[1, 1].scatter([0.0], [0.0], s=32, c="lime", edgecolor="black", lw=0.4)
    axes[1, 1].plot([-np.pi, np.pi], [np.pi, -np.pi], color="red", lw=0.8, alpha=0.7)
    axes[1, 1].set_title("Upper band: single particle reservoir")
    axes[1, 1].set_xlabel("k1")
    axes[1, 1].set_ylabel("k2")
    fig.colorbar(im1, ax=axes[1, 1], label="E - E_D")

    breakdown = particle_breakdown_metrics(params, nk=241)
    fig.suptitle(
        "Physical two-band reservoir: shifted C2T Dirac, clean hole side, particle-only reservoir "
        f"(k_D={params.root:.2f}, neck~{breakdown.neck_energy:.3f})",
        fontsize=12,
    )
    fig.savefig(out_path, dpi=230)
    plt.close(fig)
    return out_path


def plot_quick_lfd(stamp: str) -> tuple[Path, Path]:
    params = construct_single_reservoir_candidate()
    model = tb_model(params)
    q_values = np.arange(44, 145, 8)
    energy_grid = np.linspace(-0.13, 0.13, 1401)
    n_grid = np.linspace(-0.35, 0.75, 1401)

    rows = []
    b_values = []
    for q in q_values:
        spectrum = spectrum_at_flux(model, p=1, q=int(q), kmesh=(1, 7), method="dense_full")
        values = np.asarray(spectrum.eigenvalues, dtype=float).reshape(-1)
        weights = np.full(values.size, float(spectrum.state_weight_per_original_cell))
        rows.append(histogram_gaussian_dos(energy_grid, values, eta=0.0035, weights=weights))
        b_values.append(magnetic_field_tesla_from_q(int(q), p=1))
        print(f"q={int(q)} B={b_values[-1]:.4f} states={values.size}", flush=True)

    b_values = np.asarray(b_values)
    order = np.argsort(b_values)
    dos_eb = np.vstack(rows)[order]
    b_sorted = b_values[order]
    fan = regrid_dos_to_density(energy_grid, b_sorted, dos_eb, n_grid, reference_energy=0.0)

    out_dir = ROOT / "outputs" / "d3_physical_reservoir_lfd"
    fig_dir = ROOT / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(exist_ok=True)
    data_path = out_dir / f"d3_physical_reservoir_quick_lfd_{stamp}.npz"
    fig_path = fig_dir / f"d3_physical_reservoir_quick_lfd_{stamp}.png"
    np.savez_compressed(
        data_path,
        q_values=q_values,
        b_values=b_sorted,
        energy_grid=energy_grid,
        dos_EB=dos_eb,
        n_grid=n_grid,
        dos_nB=fan.dos_nB,
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.3), constrained_layout=True)
    vmax_e = np.percentile(dos_eb, 99.0)
    axes[0].pcolormesh(energy_grid, b_sorted, dos_eb, shading="auto", cmap="magma", vmin=0.0, vmax=vmax_e)
    axes[0].axvline(0.0, color="white", lw=0.8)
    axes[0].axvline(params.reservoir_energy, color="cyan", lw=0.8, ls=":")
    axes[0].set_title("Pure Peierls D(E,B)")
    axes[0].set_xlabel("E - E_D")
    axes[0].set_ylabel("B (T)")

    vmax_n = np.percentile(fan.dos_nB, 99.0)
    axes[1].pcolormesh(n_grid, b_sorted, fan.dos_nB, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax_n)
    axes[1].axvline(0.0, color="white", lw=0.8)
    axes[1].set_title("Total LFD proxy D(n,B)")
    axes[1].set_xlabel("relative filling n_rel")
    axes[1].set_ylabel("B (T)")
    fig.suptitle("Quick LFD check for shifted-C2T physical reservoir model", fontsize=12)
    fig.savefig(fig_path, dpi=230)
    plt.close(fig)
    return data_path, fig_path


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    zero_path = plot_zero_field(stamp)
    data_path, lfd_path = plot_quick_lfd(stamp)
    params = construct_single_reservoir_candidate()
    summary = {
        "params": params.__dict__,
        "reservoir_metrics": reservoir_window_metrics(params, nk=241).__dict__,
        "breakdown_metrics": particle_breakdown_metrics(params, nk=241).__dict__,
        "zero_field_figure": str(zero_path.resolve()),
        "lfd_data": str(data_path.resolve()),
        "lfd_figure": str(lfd_path.resolve()),
    }
    summary_path = ROOT / "outputs" / "d3_physical_reservoir_lfd" / f"d3_physical_reservoir_summary_{stamp}.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(zero_path.resolve())
    print(data_path.resolve())
    print(lfd_path.resolve())
    print(summary_path.resolve())


if __name__ == "__main__":
    main()
