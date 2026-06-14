"""Plot a Peierls-Hofstadter Landau-fan diagram for the D3 two-band model."""

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
    relative_filling_grid,
    slope_reversal_candidate_parameters,
    strong_bending_candidate_parameters,
    tb_model,
)
from twse2_continuum_dos.dos import (
    apply_beta_hybridization_gap,
    beta_reservoir_dos,
    magnetic_breakdown_replicated_spectrum,
    saddle_window_gaussian_dos,
    track_low_dos_branches_dp,
    regrid_dos_to_density,
)
from twse2_continuum_dos.peierls_hofstadter import (
    magnetic_field_tesla_from_q,
    rational_flux_grid_for_b_field,
    spectrum_at_flux,
)


def magnetic_dos_rows(
    params: D3TwoBandParameters,
    p_values: np.ndarray,
    q_values: np.ndarray,
    energy_grid: np.ndarray,
    *,
    saddle_energy: float,
    gamma_alpha: float,
    gamma_beta: float,
    saddle_width: float,
    beta_amplitude: float = 0.0,
    beta_center_energy: float | None = None,
    beta_width: float = 0.08,
    beta_b_star: float = 5.5,
    beta_power: float = 2.0,
    beta_onset_field: float | None = None,
    beta_onset_width: float = 0.45,
    beta_lower_energy: float | None = None,
    beta_lower_gate_width: float = 0.004,
    hybrid_gap_depth: float = 0.0,
    hybrid_gap_width: float = 0.030,
    hybrid_gap_center_slope: float = 0.0,
    hybrid_gap_onset_field: float = 4.3,
    hybrid_gap_onset_width: float = 0.35,
    reference_energy: float | None = None,
    mb_t_max: float = 0.0,
    mb_onset_field: float = 4.5,
    mb_onset_width: float = 0.5,
    mb_particle_gate_width: float = 0.01,
    kmesh: tuple[int, int] = (1, 5),
) -> tuple[np.ndarray, np.ndarray]:
    """Return B values and DOS rows from Hofstadter spectra at p/q."""
    model = tb_model(params)
    if len(p_values) != len(q_values):
        raise ValueError("p_values and q_values must have the same length")
    b_values = np.empty(len(q_values), dtype=float)
    dos_rows = np.empty((len(q_values), len(energy_grid)), dtype=float)
    for index, (p, q) in enumerate(zip(p_values, q_values)):
        spectrum = spectrum_at_flux(
            model,
            p=int(p),
            q=int(q),
            kmesh=kmesh,
            method="dense_full",
            return_eigenvectors=False,
        )
        b_values[index] = magnetic_field_tesla_from_q(int(q), p=int(p))
        values = spectrum.eigenvalues.reshape(-1)
        weights = np.full(values.size, spectrum.state_weight_per_original_cell)
        if mb_t_max > 0.0:
            if reference_energy is None:
                raise ValueError("reference_energy is required when mb_t_max is positive")
            values, weights = magnetic_breakdown_replicated_spectrum(
                values,
                weights,
                b_field=b_values[index],
                reference_energy=float(reference_energy),
                t_max=mb_t_max,
                onset_field=mb_onset_field,
                onset_width=mb_onset_width,
                particle_gate_width=mb_particle_gate_width,
                internal_degeneracy=float(params.internal_degeneracy),
                copy_degeneracy=int(params.copy_degeneracy),
            )
            degeneracy = 1.0
        else:
            degeneracy = float(params.internal_degeneracy)
        dos_rows[index] = saddle_window_gaussian_dos(
            energy_grid,
            values,
            saddle_energy=saddle_energy,
            gamma_alpha=gamma_alpha,
            gamma_beta=gamma_beta,
            saddle_width=saddle_width,
            weights=weights,
            degeneracy=degeneracy,
        )
        print(
            f"p/q={int(p):2d}/{int(q):3d} "
            f"B={b_values[index]:6.3f}T states={spectrum.eigenvalues.size:5d}"
        )
    order = np.argsort(b_values)
    b_sorted = b_values[order]
    dos_sorted = dos_rows[order]
    if beta_amplitude > 0.0:
        center = saddle_energy if beta_center_energy is None else float(beta_center_energy)
        dos_sorted = dos_sorted + beta_reservoir_dos(
            energy_grid,
            b_sorted,
            center_energy=center,
            width=beta_width,
            amplitude=beta_amplitude,
            b_star=beta_b_star,
            power=beta_power,
            onset_field=beta_onset_field,
            onset_width=beta_onset_width,
            lower_energy=beta_lower_energy,
            lower_gate_width=beta_lower_gate_width,
        )
    if hybrid_gap_depth > 0.0:
        center = saddle_energy if beta_center_energy is None else float(beta_center_energy)
        dos_sorted = apply_beta_hybridization_gap(
            energy_grid,
            b_sorted,
            dos_sorted,
            center_energy=center,
            width=hybrid_gap_width,
            depth=hybrid_gap_depth,
            onset_field=hybrid_gap_onset_field,
            onset_width=hybrid_gap_onset_width,
            center_slope=hybrid_gap_center_slope,
            center_reference_field=hybrid_gap_onset_field,
        )
    return b_sorted, dos_sorted


def selected_parameters(name: str) -> D3TwoBandParameters:
    if name == "strong":
        return strong_bending_candidate_parameters()
    if name == "slope_reversal":
        return slope_reversal_candidate_parameters()
    raise ValueError(f"unknown candidate: {name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", choices=["strong", "slope_reversal"], default="slope_reversal")
    parser.add_argument("--mode", choices=["mb_spacing", "bending"], default="mb_spacing")
    parser.add_argument("--b-min", type=float, default=0.76)
    parser.add_argument("--b-max", type=float, default=10.0)
    parser.add_argument("--b-points", type=int, default=120)
    parser.add_argument("--q-max", type=int, default=263)
    parser.add_argument("--show-branches", action="store_true")
    args = parser.parse_args()

    params = selected_parameters(args.candidate)
    report = candidate_report(params, nk=151)
    e_dirac = float(report["selected_dirac"]["energy"])

    p_values, q_values, sampled_b_values = rational_flux_grid_for_b_field(
        b_min=args.b_min,
        b_max=args.b_max,
        n_points=args.b_points,
        q_max=args.q_max,
    )
    energy_grid = np.linspace(e_dirac - 1.0, e_dirac + 1.0, 4001)
    n_grid = relative_filling_grid(-1.0, 1.0, 2401)
    saddle_energy = float(report["upper_saddles_above_dirac"][0]["energy"])
    gamma_alpha = 0.003
    if args.mode == "mb_spacing":
        gamma_alpha = 0.0035
        gamma_beta = 0.0035
        saddle_width = 0.030
        beta_amplitude = 0.0
        hybrid_gap_depth = 0.0
        mb_t_max = 0.055
        mb_onset_field = 4.5
        mb_onset_width = 0.55
        mb_particle_gate_width = 0.008
    else:
        gamma_alpha = 0.003
        gamma_beta = 0.045
        saddle_width = 0.050
        beta_amplitude = 1.10
        hybrid_gap_depth = 0.65
        mb_t_max = 0.0
        mb_onset_field = 4.5
        mb_onset_width = 0.55
        mb_particle_gate_width = 0.008
    beta_center_energy = saddle_energy
    beta_width = 0.085
    beta_b_star = 4.5
    beta_power = 2.0
    beta_onset_field = 4.2
    beta_onset_width = 0.32
    beta_lower_energy = e_dirac + 0.005
    beta_lower_gate_width = 0.004
    hybrid_gap_width = 0.030
    hybrid_gap_center_slope = 0.015
    hybrid_gap_onset_field = 4.3
    hybrid_gap_onset_width = 0.35
    show_dark_branch_overlay = bool(args.show_branches)

    b_values, dos_eb = magnetic_dos_rows(
        params,
        p_values,
        q_values,
        energy_grid,
        saddle_energy=saddle_energy,
        gamma_alpha=gamma_alpha,
        gamma_beta=gamma_beta,
        saddle_width=saddle_width,
        beta_amplitude=beta_amplitude,
        beta_center_energy=beta_center_energy,
        beta_width=beta_width,
        beta_b_star=beta_b_star,
        beta_power=beta_power,
        beta_onset_field=beta_onset_field,
        beta_onset_width=beta_onset_width,
        beta_lower_energy=beta_lower_energy,
        beta_lower_gate_width=beta_lower_gate_width,
        hybrid_gap_depth=hybrid_gap_depth,
        hybrid_gap_width=hybrid_gap_width,
        hybrid_gap_center_slope=hybrid_gap_center_slope,
        hybrid_gap_onset_field=hybrid_gap_onset_field,
        hybrid_gap_onset_width=hybrid_gap_onset_width,
        reference_energy=e_dirac,
        mb_t_max=mb_t_max,
        mb_onset_field=mb_onset_field,
        mb_onset_width=mb_onset_width,
        mb_particle_gate_width=mb_particle_gate_width,
        kmesh=(1, 13),
    )
    fan = regrid_dos_to_density(energy_grid, b_values, dos_eb, n_grid, reference_energy=e_dirac)

    out_dir = ROOT / "figures"
    out_dir.mkdir(exist_ok=True)
    data_dir = ROOT / "outputs" / "d3_two_band_lfd"
    data_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out_path = out_dir / f"d3_two_band_peierls_lfd_{stamp}.png"
    data_path = data_dir / f"d3_two_band_peierls_lfd_{stamp}.npz"

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

    dos_n_vmax = max(np.percentile(fan.dos_nB, 99.0), 1e-12)
    mesh_n = axes[1].pcolormesh(
        n_grid,
        b_values,
        fan.dos_nB,
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=dos_n_vmax,
    )
    if show_dark_branch_overlay:
        dark_branches = track_low_dos_branches_dp(
            n_grid,
            b_values,
            fan.dos_nB,
            seed_n_values=[0.08, 0.13, 0.18, 0.24, 0.31],
            seed_slopes=[0.040, 0.055, 0.070, 0.085, 0.100],
            n_min=0.02,
            n_max=0.98,
            percentile_cut=60.0,
            candidates_per_row=16,
            seed_weight=300.0,
            continuity_weight=70.0,
            slope_weight=10.0,
            max_step=0.10,
        )
        for branch in dark_branches:
            finite = np.isfinite(branch["n"])
            axes[1].plot(branch["n"][finite], branch["b"][finite], color="white", lw=1.05, alpha=0.86)
    axes[1].axvline(0.0, color="white", lw=0.9, alpha=0.85)
    axes[1].set_xlim(-1.0, 1.0)
    axes[1].set_title("Density/filling fan: D(n_rel, B)")
    axes[1].set_xlabel("relative filling n_rel, Dirac half filling at 0")
    axes[1].set_ylabel("B (T)")
    fig.colorbar(mesh_n, ax=axes[1], label="DOS at mu(nu,B)")

    fig.suptitle(
        f"D3 Peierls-Hofstadter LFD ({args.candidate}, {args.mode})\n"
        f"B={sampled_b_values.min():.2f}..{sampled_b_values.max():.2f}T, "
        f"flux points={len(q_values)}, q<= {q_values.max()}, kmesh=(1,13), "
        rf"$\Gamma_\alpha={gamma_alpha:g},\Gamma_\beta={gamma_beta:g},w_s={saddle_width:g},A_\beta={beta_amplitude:g}$",
        fontsize=11,
    )
    fig.savefig(out_path, dpi=260, bbox_inches="tight")
    plt.close(fig)

    np.savez_compressed(
        data_path,
        candidate=np.asarray(args.candidate),
        p_values=p_values,
        q_values=q_values,
        sampled_b_values=sampled_b_values,
        energy_grid=energy_grid,
        b_values=b_values,
        dos_eb=dos_eb,
        n_grid=n_grid,
        dos_nb=fan.dos_nB,
        mu_nb=fan.mu_nB,
        n_eb=fan.n_EB,
        e_dirac=np.asarray(e_dirac),
        saddle_energy=np.asarray(saddle_energy),
        gamma_alpha=np.asarray(gamma_alpha),
        gamma_beta=np.asarray(gamma_beta),
        saddle_width=np.asarray(saddle_width),
        beta_amplitude=np.asarray(beta_amplitude),
        beta_center_energy=np.asarray(beta_center_energy),
        beta_width=np.asarray(beta_width),
        beta_onset_field=np.asarray(beta_onset_field),
        beta_onset_width=np.asarray(beta_onset_width),
        hybrid_gap_depth=np.asarray(hybrid_gap_depth),
        hybrid_gap_width=np.asarray(hybrid_gap_width),
        hybrid_gap_center_slope=np.asarray(hybrid_gap_center_slope),
        hybrid_gap_onset_field=np.asarray(hybrid_gap_onset_field),
        hybrid_gap_onset_width=np.asarray(hybrid_gap_onset_width),
        mode=np.asarray(args.mode),
        mb_t_max=np.asarray(mb_t_max),
        mb_onset_field=np.asarray(mb_onset_field),
        mb_onset_width=np.asarray(mb_onset_width),
        mb_particle_gate_width=np.asarray(mb_particle_gate_width),
    )

    print(out_path.resolve())
    print(data_path.resolve())
    print(f"candidate={args.candidate}")
    print(f"mode={args.mode}")
    print(f"flux_count={len(q_values)}")
    print(f"q_range=({q_values.min()}, {q_values.max()})")
    print(f"target_B_range=({args.b_min:.6f}, {args.b_max:.6f})")
    print(f"E_D={e_dirac:.9f}")
    print(f"nearest_upper_saddle_dE={report['upper_saddles_above_dirac'][0]['energy_minus_dirac']:.9f}")
    print(
        "beta_reservoir="
        f"(A={beta_amplitude:.6f}, center_minus_ED={beta_center_energy - e_dirac:.6f}, "
        f"width={beta_width:.6f}, B_star={beta_b_star:.6f}, power={beta_power:.6f}, "
        f"B_on={beta_onset_field:.6f}, B_on_width={beta_onset_width:.6f})"
    )
    print(
        "hybrid_gap="
        f"(depth={hybrid_gap_depth:.6f}, width={hybrid_gap_width:.6f}, "
        f"center_slope={hybrid_gap_center_slope:.6f}, B_on={hybrid_gap_onset_field:.6f})"
    )
    print(
        "magnetic_breakdown="
        f"(t_max={mb_t_max:.6f}, B_on={mb_onset_field:.6f}, "
        f"B_on_width={mb_onset_width:.6f}, particle_gate_width={mb_particle_gate_width:.6f})"
    )
    print("density_axis=relative filling n_rel with Dirac half filling at 0")
    print(f"B_range=({b_values.min():.6f}, {b_values.max():.6f})")
    b_steps = np.diff(b_values)
    print(
        "B_step_stats="
        f"(min={b_steps.min():.6f}, median={np.median(b_steps):.6f}, "
        f"max={b_steps.max():.6f}, std={b_steps.std():.6f})"
    )
    print(f"energy_window=({energy_grid[0] - e_dirac:.6f}, {energy_grid[-1] - e_dirac:.6f}) relative to E_D")
    print(f"n_rel_window=({n_grid[0]:.6f}, {n_grid[-1]:.6f})")
    print(f"dark_branch_overlay={'dp' if show_dark_branch_overlay else 'off'}")


if __name__ == "__main__":
    main()
