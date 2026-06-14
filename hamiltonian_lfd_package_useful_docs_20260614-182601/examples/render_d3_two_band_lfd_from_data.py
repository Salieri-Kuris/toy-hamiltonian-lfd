"""Render a D3 two-band Landau-fan figure from a saved NPZ data file."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _scalar_text(data: np.lib.npyio.NpzFile, key: str, default: str) -> str:
    if key not in data:
        return default
    value = data[key]
    if value.shape == ():
        return str(value.item())
    return str(value)


def render_lfd(data_path: Path, output_path: Path | None = None) -> Path:
    with np.load(data_path, allow_pickle=False) as data:
        energy_grid = np.asarray(data["energy_grid"], dtype=float)
        b_values = np.asarray(data["b_values"], dtype=float)
        dos_eb = np.asarray(data["dos_eb"], dtype=float)
        n_grid = np.asarray(data["n_grid"], dtype=float)
        dos_nb = np.asarray(data["dos_nb"], dtype=float)
        e_dirac = float(np.asarray(data["e_dirac"]).item())
        candidate = _scalar_text(data, "candidate", "unknown")
        mode = _scalar_text(data, "mode", "unknown")
        q_values = np.asarray(data["q_values"], dtype=int) if "q_values" in data else np.asarray([], dtype=int)
        gamma_alpha = float(np.asarray(data["gamma_alpha"]).item()) if "gamma_alpha" in data else np.nan
        gamma_beta = float(np.asarray(data["gamma_beta"]).item()) if "gamma_beta" in data else np.nan
        saddle_width = float(np.asarray(data["saddle_width"]).item()) if "saddle_width" in data else np.nan
        beta_amplitude = float(np.asarray(data["beta_amplitude"]).item()) if "beta_amplitude" in data else np.nan

    if output_path is None:
        out_dir = ROOT / "figures"
        out_dir.mkdir(exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        output_path = out_dir / f"{data_path.stem}_render_{stamp}.png"
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.9), constrained_layout=True)

    dos_e_vmax = max(np.percentile(dos_eb, 99.0), 1e-12)
    mesh_e = axes[0].pcolormesh(
        energy_grid - e_dirac,
        b_values,
        dos_eb,
        shading="auto",
        cmap="magma",
        vmin=0.0,
        vmax=dos_e_vmax,
    )
    axes[0].axvline(0.0, color="white", lw=0.8, alpha=0.75)
    axes[0].set_xlim(-1.0, 1.0)
    axes[0].set_title("Peierls-Hofstadter DOS: D(E, B)")
    axes[0].set_xlabel("E - E_D")
    axes[0].set_ylabel("B (T)")
    fig.colorbar(mesh_e, ax=axes[0], label="DOS")

    dos_n_vmax = max(np.percentile(dos_nb, 99.0), 1e-12)
    mesh_n = axes[1].pcolormesh(
        n_grid,
        b_values,
        dos_nb,
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=dos_n_vmax,
    )
    axes[1].axvline(0.0, color="white", lw=0.9, alpha=0.85)
    axes[1].set_xlim(-1.0, 1.0)
    axes[1].set_title("Density/filling fan: D(n_rel, B)")
    axes[1].set_xlabel("relative filling n_rel, Dirac half filling at 0")
    axes[1].set_ylabel("B (T)")
    fig.colorbar(mesh_n, ax=axes[1], label="DOS at mu(nu,B)")

    q_text = f", q<= {q_values.max()}" if q_values.size else ""
    fig.suptitle(
        f"D3 Peierls-Hofstadter LFD ({candidate}, {mode})\n"
        f"B={b_values.min():.2f}..{b_values.max():.2f}T, "
        f"flux points={len(b_values)}{q_text}, "
        rf"$\Gamma_\alpha={gamma_alpha:g},\Gamma_\beta={gamma_beta:g},w_s={saddle_width:g},A_\beta={beta_amplitude:g}$",
        fontsize=11,
    )
    fig.savefig(output_path, dpi=260, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="Saved NPZ file from plot_d3_two_band_lfd.py")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    out_path = render_lfd(args.data, args.output)
    print(out_path.resolve())
    print(f"data_source={args.data.resolve()}")


if __name__ == "__main__":
    main()
