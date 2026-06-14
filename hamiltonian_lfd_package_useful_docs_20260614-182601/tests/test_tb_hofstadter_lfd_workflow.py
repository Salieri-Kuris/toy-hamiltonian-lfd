import numpy as np

from twse2_continuum_dos.tb_hofstadter_lfd_workflow import (
    compute_pure_hofstadter_dos,
    selected_d3_parameters,
)


def test_selected_d3_parameters_accepts_named_manual_candidates():
    assert selected_d3_parameters("strong") is not None
    assert selected_d3_parameters("slope_reversal") is not None


def test_compute_pure_hofstadter_dos_returns_sorted_total_dos_and_fan():
    params = selected_d3_parameters("strong")
    energy_grid = np.linspace(-0.25, 0.25, 151)
    n_grid = np.linspace(-0.2, 0.2, 101)

    result = compute_pure_hofstadter_dos(
        params,
        q_values=np.array([9, 7]),
        energy_grid=energy_grid,
        n_grid=n_grid,
        eta=0.01,
        kmesh=(1, 2),
        reference_energy=0.0,
    )

    assert result.b_values.shape == (2,)
    assert np.all(np.diff(result.b_values) > 0.0)
    assert result.dos_EB.shape == (2, len(energy_grid))
    assert result.fan.dos_nB.shape == (2, len(n_grid))
    assert result.q_values.tolist() == [9, 7]
    assert np.all(result.dos_EB >= 0.0)
