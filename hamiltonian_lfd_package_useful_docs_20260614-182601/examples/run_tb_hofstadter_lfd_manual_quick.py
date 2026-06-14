"""Quick pure-Peierls TB Hofstadter LFD run following the manual checklist."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from twse2_continuum_dos.d3_two_band import candidate_report, relative_filling_grid
from twse2_continuum_dos.d3_two_band import D3TwoBandParameters
from twse2_continuum_dos.tb_hofstadter_lfd_workflow import (
    compute_pure_hofstadter_dos,
    selected_d3_parameters,
)


def q_grid(q_min: int, q_max: int, q_step: int) -> np.ndarray:
    """Return denominators for a quick p=1 Hofstadter scan."""
    if q_min <= 0 or q_max <= 0 or q_step <= 0:
        raise ValueError("q_min, q_max, and q_step must be positive")
    if q_min > q_max:
        raise ValueError("q_min must not exceed q_max")
    return np.arange(int(q_min), int(q_max) + 1, int(q_step), dtype=int)


def to_jsonable_report(report: dict) -> dict:
    """Keep the diagnostic report compact and JSON serializable."""
    return {
        "parameters": report["parameters"],
        "selected_dirac": report["selected_dirac"],
        "nu_offset": report["nu_offset"],
        "selected_dirac_filling": report["selected_dirac_filling"],
        "upper_saddles_above_dirac": report["upper_saddles_above_dirac"][:4],
        "lower_saddles_below_dirac": report["lower_saddles_below_dirac"][:4],
    }


def parameters_from_args(args: argparse.Namespace) -> D3TwoBandParameters:
    """Return selected D3 parameters, optionally overridden from the CLI."""
    params = selected_d3_parameters(args.candidate)
    m_values = [params.m0, params.m1, params.m2, params.m3] if args.m is None else args.m
    if len(m_values) != 4:
        raise ValueError("--m requires exactly four values for m0..m3")
    lam_value = params.lam if args.lam is None else float(args.lam)
    u_values = [params.u1, params.u2, params.u3, params.u4, params.u5, params.u6, params.u7] if args.u is None else args.u
    if len(u_values) != 7:
        raise ValueError("--u requires exactly seven values for u1..u7")
    return D3TwoBandParameters(
        m0=m_values[0],
        m1=m_values[1],
        m2=m_values[2],
        m3=m_values[3],
        lam=lam_value,
        u0=params.u0,
        u1=u_values[0],
        u2=u_values[1],
        u3=u_values[2],
        u4=u_values[3],
        u5=u_values[4],
        u6=u_values[5],
        u7=u_values[6],
        internal_degeneracy=params.internal_degeneracy,
        copy_degeneracy=params.copy_degeneracy,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", choices=["strong", "slope_reversal"], default="slope_reversal")
    parser.add_argument("--m", type=float, nargs="*", help="optional m0..m3 override for changed mass texture")
    parser.add_argument("--lam", type=float, help="optional lambda override for changed off-diagonal strength")
    parser.add_argument("--u", type=float, nargs="*", help="optional u1..u7 override for changed d0 Hamiltonian")
    parser.add_argument("--q-min", type=int, default=24)
    parser.add_argument("--q-max", type=int, default=120)
    parser.add_argument("--q-step", type=int, default=6)
    parser.add_argument("--kmesh2", type=int, default=5)
    parser.add_argument("--eta", type=float, default=0.006)
    parser.add_argument("--energy-half-width", type=float, default=0.18)
    parser.add_argument("--energy-points", type=int, default=1801)
    parser.add_argument("--n-min", type=float, default=-0.6)
    parser.add_argument("--n-max", type=float, default=1.0)
    parser.add_argument("--n-points", type=int, default=1601)
    args = parser.parse_args()

    params = parameters_from_args(args)
    report = candidate_report(params, nk=121)
    e_dirac = float(report["selected_dirac"]["energy"])
    saddle_energy = float(report["upper_saddles_above_dirac"][0]["energy"])

    energy_grid = np.linspace(
        e_dirac - float(args.energy_half_width),
        e_dirac + float(args.energy_half_width),
        int(args.energy_points),
    )
    n_grid = relative_filling_grid(float(args.n_min), float(args.n_max), int(args.n_points))
    q_values = q_grid(args.q_min, args.q_max, args.q_step)
    result = compute_pure_hofstadter_dos(
        params,
        q_values=q_values,
        energy_grid=energy_grid,
        n_grid=n_grid,
        eta=float(args.eta),
        kmesh=(1, int(args.kmesh2)),
        reference_energy=e_dirac,
    )

    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    candidate_label = args.candidate if args.m is None and args.lam is None and args.u is None else f"{args.candidate}_override"
    figure_dir = ROOT / "figures"
    data_dir = ROOT / "outputs" / "tb_hofstadter_lfd_manual_quick"
    figure_dir.mkdir(exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_path = figure_dir / f"tb_hofstadter_lfd_manual_quick_{candidate_label}_{stamp}.png"
    data_path = data_dir / f"tb_hofstadter_lfd_manual_quick_{candidate_label}_{stamp}.npz"
    summary_path = data_dir / f"tb_hofstadter_lfd_manual_quick_{candidate_label}_{stamp}.json"

    np.savez_compressed(
        data_path,
        candidate=np.asarray(candidate_label),
        q_values=result.q_values,
        b_values=result.b_values,
        energy_grid=result.energy_grid,
        dos_EB=result.dos_EB,
        n_grid=result.fan.n_grid,
        dos_nB=result.fan.dos_nB,
        mu_nB=result.fan.mu_nB,
        n_EB=result.fan.n_EB,
        eta=np.asarray(result.eta),
        kmesh=np.asarray(result.kmesh, dtype=int),
        e_dirac=np.asarray(e_dirac),
        saddle_energy=np.asarray(saddle_energy),
    )

    summary = {
        "candidate": candidate_label,
        "base_candidate": args.candidate,
        "parameters": params.__dict__,
        "figure_path": str(figure_path.resolve()),
        "data_path": str(data_path.resolve()),
        "q_count": int(len(result.q_values)),
        "q_range": [int(result.q_values.min()), int(result.q_values.max())],
        "b_range_tesla": [float(result.b_values.min()), float(result.b_values.max())],
        "eta": float(result.eta),
        "kmesh": list(result.kmesh),
        "e_dirac": e_dirac,
        "nearest_upper_saddle_minus_dirac": float(saddle_energy - e_dirac),
        "manual_constraint": "pure Peierls total DOS; no beta_reservoir_dos, no hand-drawn fan lines",
        "zero_field_report": to_jsonable_report(report),
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), constrained_layout=True)
    vmax_e = max(float(np.percentile(result.dos_EB, 99.2)), 1e-12)
    mesh_e = axes[0].pcolormesh(
        result.energy_grid - e_dirac,
        result.b_values,
        result.dos_EB,
        shading="auto",
        cmap="magma",
        vmin=0.0,
        vmax=vmax_e,
    )
    axes[0].axvline(0.0, color="white", lw=0.8, alpha=0.75)
    axes[0].axvline(saddle_energy - e_dirac, color="cyan", lw=0.9, alpha=0.9)
    axes[0].set_xlabel("E - E_D")
    axes[0].set_ylabel("B (T)")
    axes[0].set_title("Pure Peierls total DOS: D(E, B)")
    fig.colorbar(mesh_e, ax=axes[0], label="DOS")

    vmax_n = max(float(np.percentile(result.fan.dos_nB, 99.2)), 1e-12)
    mesh_n = axes[1].pcolormesh(
        result.fan.n_grid,
        result.b_values,
        result.fan.dos_nB,
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax_n,
    )
    axes[1].axvline(0.0, color="white", lw=0.8, alpha=0.75)
    axes[1].set_xlabel("relative filling n_rel, E_D at 0")
    axes[1].set_ylabel("B (T)")
    axes[1].set_title("Total LFD proxy: D(n, B)")
    fig.colorbar(mesh_n, ax=axes[1], label="DOS at mu(n,B)")

    fig.suptitle(
        f"TB manual quick check: {candidate_label}, "
        f"q={result.q_values.min()}..{result.q_values.max()}, "
        f"kmesh={result.kmesh}, eta={result.eta:g}, "
        f"upper saddle - E_D={saddle_energy - e_dirac:.5f}",
        fontsize=11,
    )
    fig.savefig(figure_path, dpi=240, bbox_inches="tight")
    plt.close(fig)

    print(figure_path.resolve())
    print(data_path.resolve())
    print(summary_path.resolve())
    print(f"candidate={candidate_label}")
    print(f"q_count={len(result.q_values)}")
    print(f"B_range=({result.b_values.min():.6f}, {result.b_values.max():.6f}) T")
    print(f"E_D={e_dirac:.9f}")
    print(f"nearest_upper_saddle_dE={saddle_energy - e_dirac:.9f}")
    print("manual_constraint=pure Peierls total DOS; no manual beta reservoir")


if __name__ == "__main__":
    main()
