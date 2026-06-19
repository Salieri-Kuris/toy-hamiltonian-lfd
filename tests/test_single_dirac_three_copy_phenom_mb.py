import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from plot_single_dirac_three_copy_phenom_mb import (
    PhenomParams,
    build_dos_data,
    mb_diagnostics,
    split_three_copy_levels,
)


def test_default_parameters_make_particle_mb_stronger_than_hole_mb():
    params = PhenomParams()

    diagnostics = mb_diagnostics(params, energy=params.Ebeta, B=0.75)

    assert diagnostics["P_particle"] > 0.10
    assert diagnostics["P_hole"] < 1.0e-2
    assert diagnostics["particle_split"] > params.Gamma_LL
    assert diagnostics["hole_split"] < params.Gamma_LL


def test_zeroth_landau_level_is_present_at_dirac_energy_without_splitting():
    params = PhenomParams(alpha=0.0, t_scale=0.0)

    particle_levels, _weights = split_three_copy_levels(0.7, "particle", params)
    hole_levels, _weights = split_three_copy_levels(0.7, "hole", params)

    assert np.count_nonzero(np.isclose(particle_levels, params.ED)) == 3
    assert np.count_nonzero(np.isclose(hole_levels, params.ED)) == 3


def test_build_dos_data_returns_finite_nonnegative_grids():
    params = PhenomParams(Nmax=12)
    B_values = np.linspace(0.12, 0.80, 9)
    energy_grid = np.linspace(-1.2, 1.2, 401)
    n_grid = np.linspace(-0.18, 0.18, 241)

    data = build_dos_data(
        params,
        B_values=B_values,
        energy_grid=energy_grid,
        n_grid=n_grid,
    )

    assert data.source_dos_EB.shape == (len(B_values), len(energy_grid))
    assert data.fan.dos_nB.shape == (len(B_values), len(n_grid))
    assert np.all(np.isfinite(data.source_dos_EB))
    assert np.nanmin(data.source_dos_EB) >= 0.0
    assert np.nanmin(data.fan.dos_nB) >= 0.0
