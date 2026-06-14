"""Compute and cache a higher-resolution D3 two-band Hofstadter LFD dataset."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from twse2_continuum_dos.d3_two_band import strong_bending_candidate_parameters, candidate_report, relative_filling_grid, tb_model
from twse2_continuum_dos.dos import saddle_window_gaussian_dos, regrid_dos_to_density
from twse2_continuum_dos.peierls_hofstadter import magnetic_field_tesla_from_q, spectrum_at_flux


def q_grid() -> np.ndarray:
    """Return an ultra-dense rational-flux grid sized for a sub-hour local run."""
    low = np.arange(20, 301, 1)
    mid = np.arange(304, 601, 8)
    high = np.arange(608, 801, 16)
    return np.unique(np.concatenate((low, mid, high))).astype(int)


def compute_spectra(q_values: np.ndarray, kmesh: tuple[int, int]) -> dict[str, np.ndarray]:
    """Compute Hofstadter spectra and pack ragged eigenvalue arrays."""
    params = strong_bending_candidate_parameters()
    model = tb_model(params)
    report = candidate_report(params, nk=151)
    e_dirac = float(report["selected_dirac"]["energy"])
    saddle = report["upper_saddles_above_dirac"][0]
    saddle_energy = float(saddle["energy"])
    saddle_de = float(saddle["energy_minus_dirac"])

    all_values: list[np.ndarray] = []
    offsets = [0]
    b_values = []
    state_weights = []
    for q in q_values:
        spectrum = spectrum_at_flux(
            model,
            p=1,
            q=int(q),
            kmesh=kmesh,
            method="dense_full",
            return_eigenvectors=False,
        )
        values = np.asarray(spectrum.eigenvalues, dtype=float).reshape(-1)
        all_values.append(values)
        offsets.append(offsets[-1] + len(values))
        b_values.append(magnetic_field_tesla_from_q(int(q), p=1))
        state_weights.append(float(spectrum.state_weight_per_original_cell))
        print(
            f"q={int(q):3d} B={b_values[-1]:7.4f}T "
            f"dim={2 * int(q):4d} kpts={len(spectrum.k_points):2d} states={len(values):6d}",
            flush=True,
        )

    return {
        "q_values": q_values.astype(int),
        "b_values": np.asarray(b_values, dtype=float),
        "eigenvalues_flat": np.concatenate(all_values),
        "offsets": np.asarray(offsets, dtype=int),
        "state_weight_per_original_cell": np.asarray(state_weights, dtype=float),
        "e_dirac": np.asarray(e_dirac),
        "nearest_upper_saddle_energy": np.asarray(saddle_energy),
        "nearest_upper_saddle_de": np.asarray(saddle_de),
        "kmesh": np.asarray(kmesh, dtype=int),
    }


def row_values(data: dict[str, np.ndarray], index: int) -> np.ndarray:
    start = int(data["offsets"][index])
    stop = int(data["offsets"][index + 1])
    return data["eigenvalues_flat"][start:stop]


def dos_for_broadening(
    data: dict[str, np.ndarray],
    energy_grid: np.ndarray,
    *,
    gamma_alpha: float,
    gamma_beta: float,
    saddle_width: float,
) -> np.ndarray:
    params = strong_bending_candidate_parameters()
    saddle_energy = float(data["nearest_upper_saddle_energy"])
    rows = np.empty((len(data["q_values"]), len(energy_grid)), dtype=float)
    for i in range(len(data["q_values"])):
        values = row_values(data, i)
        weights = np.full(len(values), float(data["state_weight_per_original_cell"][i]))
        rows[i] = saddle_window_gaussian_dos(
            energy_grid,
            values,
            saddle_energy=saddle_energy,
            gamma_alpha=gamma_alpha,
            gamma_beta=gamma_beta,
            saddle_width=saddle_width,
            weights=weights,
            degeneracy=float(params.internal_degeneracy),
        )
    order = np.argsort(data["b_values"])
    return rows[order]


def plot_broadening_grid(data: dict[str, np.ndarray], out_dir: Path, stamp: str) -> Path:
    e_dirac = float(data["e_dirac"])
    b_values = np.asarray(data["b_values"], dtype=float)
    order = np.argsort(b_values)
    b_sorted = b_values[order]
    energy_grid = np.linspace(e_dirac - 0.24, e_dirac + 0.24, 3201)
    n_grid = relative_filling_grid(-0.55, 1.35, 2601)
    broadenings = [
        (0.003, 0.010, 0.025),
        (0.003, 0.012, 0.030),
        (0.003, 0.016, 0.040),
        (0.004, 0.018, 0.045),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.0), constrained_layout=True, sharex=True, sharey=True)
    for ax, (gamma_alpha, gamma_beta, saddle_width) in zip(axes.ravel(), broadenings):
        dos_eb = dos_for_broadening(
            data,
            energy_grid,
            gamma_alpha=gamma_alpha,
            gamma_beta=gamma_beta,
            saddle_width=saddle_width,
        )
        fan = regrid_dos_to_density(energy_grid, b_sorted, dos_eb, n_grid, reference_energy=e_dirac)
        mesh = ax.pcolormesh(
            n_grid,
            b_sorted,
            fan.dos_nB,
            shading="auto",
            cmap="viridis",
            vmin=0.0,
            vmax=max(np.percentile(fan.dos_nB, 99.0), 1e-12),
        )
        ax.axvline(0.0, color="white", lw=0.8, alpha=0.85)
        ax.set_title(rf"$\Gamma_\alpha={gamma_alpha:g},\Gamma_\beta={gamma_beta:g},w_s={saddle_width:g}$")
        ax.set_xlabel("n_rel, Dirac half filling at 0")
        ax.set_ylabel("B (T)")
        fig.colorbar(mesh, ax=ax, label="DOS")

    fig.suptitle(
        "Tuned D3 two-band particle-side LFD vs broadening "
        f"(q={data['q_values'].min()}..{data['q_values'].max()}, kmesh={tuple(data['kmesh'])})",
        fontsize=12,
    )
    out_path = out_dir / f"d3_two_band_highres_broadening_grid_{stamp}.png"
    fig.savefig(out_path, dpi=230)
    plt.close(fig)
    return out_path


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = ROOT / "outputs" / "d3_two_band_highres_lfd"
    figure_dir = ROOT / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(exist_ok=True)

    q_values = q_grid()
    kmesh = (1, 21)
    data = compute_spectra(q_values, kmesh)
    data_path = out_dir / f"d3_two_band_highres_spectra_{stamp}.npz"
    np.savez_compressed(data_path, **data)

    figure_path = plot_broadening_grid(data, figure_dir, stamp)

    print(data_path.resolve())
    print(figure_path.resolve())
    print(f"E_D={float(data['e_dirac']):.9f}")
    print(f"nearest_upper_saddle_dE={float(data['nearest_upper_saddle_de']):.9f}")
    print(f"q_count={len(q_values)}")
    print(f"B_range=({data['b_values'].min():.6f}, {data['b_values'].max():.6f})")


if __name__ == "__main__":
    main()
