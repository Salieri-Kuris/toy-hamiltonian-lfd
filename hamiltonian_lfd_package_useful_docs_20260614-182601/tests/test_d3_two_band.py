import numpy as np

from twse2_continuum_dos.d3_two_band import (
    D3TwoBandParameters,
    background_candidate_parameters,
    slope_reversal_candidate_parameters,
    strong_bending_candidate_parameters,
    starting_candidate_parameters,
    calibrate_nu_offset,
    candidate_report,
    filling_at_energy,
    find_line_dirac_roots,
    g_harmonic,
    hamiltonian,
    line_spectrum,
    relative_filling_grid,
    tb_model,
)


def test_g_harmonic_vanishes_on_representative_c2t_line():
    k1 = np.linspace(-np.pi, np.pi, 51)

    assert np.allclose(g_harmonic(k1, -k1), 0.0, atol=1e-14)


def test_g_harmonic_does_not_vanish_on_unitary_c2_line():
    k1 = np.array([-1.7, -0.8, 0.9, 1.6])

    assert not np.allclose(g_harmonic(k1, 0.0), 0.0, atol=1e-6)


def test_g_harmonic_vanishes_on_c2t_boundary_line():
    k1 = np.linspace(-np.pi, np.pi, 51)

    assert np.allclose(g_harmonic(k1, np.pi - k1), 0.0, atol=1e-14)


def test_d3_two_band_hamiltonian_is_hermitian():
    params = D3TwoBandParameters()

    h = hamiltonian(0.37, -1.11, params)

    assert h.shape == (2, 2)
    assert np.allclose(h, h.conj().T, atol=1e-14)


def test_default_parameters_are_current_starting_hamiltonian():
    assert D3TwoBandParameters() == starting_candidate_parameters()
    assert D3TwoBandParameters() != background_candidate_parameters()


def test_tb_model_reconstructs_direct_d3_two_band_hamiltonian():
    params = D3TwoBandParameters(u5=0.017, u6=-0.011, u7=0.006)
    model = tb_model(params)

    assert model.n_orb == 2
    for k1, k2 in [(0.13, 0.29), (1.2, -0.4), (-2.0, 2.4)]:
        assert np.allclose(model.bloch_hamiltonian(k1, k2), hamiltonian(k1, k2, params), atol=1e-12)
    assert model.hermiticity_error() < 1e-12


def test_saved_two_band_candidates_are_two_bands_per_valley():
    for params in [
        starting_candidate_parameters(),
        background_candidate_parameters(),
        strong_bending_candidate_parameters(),
        slope_reversal_candidate_parameters(),
    ]:
        model = tb_model(params)
        assert model.n_orb == 2
        assert hamiltonian(0.21, -0.37, params).shape == (2, 2)


def test_default_candidate_has_line_dirac_root_with_nonzero_velocities():
    params = D3TwoBandParameters()

    roots = find_line_dirac_roots(params)
    report = candidate_report(params, nk=81)

    assert roots
    assert np.isclose(report["selected_dirac"]["k1"] + report["selected_dirac"]["k2"], 0.0, atol=1e-9)
    assert abs(g_harmonic(report["selected_dirac"]["k1"], report["selected_dirac"]["k2"])) < 1e-12
    assert abs(g_harmonic(report["selected_dirac"]["k1"], 0.0)) > 1e-3
    assert abs(report["selected_dirac"]["v_parallel"]) > 1e-3
    assert abs(report["selected_dirac"]["v_perp"]) > 1e-3


def test_filling_offset_calibrates_selected_dirac_to_nu_minus_two():
    params = D3TwoBandParameters()
    report = candidate_report(params, nk=71)
    e_dirac = report["selected_dirac"]["energy"]

    offset = calibrate_nu_offset(params, e_dirac, target_nu=-2.0, nk=71)

    assert np.isclose(filling_at_energy(params, e_dirac, nu_offset=offset, nk=71), -2.0)


def test_default_candidate_has_particle_side_upper_saddle():
    report = candidate_report(D3TwoBandParameters(), nk=101)
    saddles = report["upper_saddles_above_dirac"]

    assert saddles
    assert saddles[0]["energy_minus_dirac"] > 0.0
    assert saddles[0]["hessian_det"] < 0.0


def test_d3_two_band_public_exports_are_available():
    import twse2_continuum_dos as tw

    params = tw.D3TwoBandParameters()
    report = tw.d3_candidate_report(params, nk=51)

    assert "selected_dirac" in report
    assert tw.d3_strong_bending_candidate_parameters() == strong_bending_candidate_parameters()
    assert tw.d3_slope_reversal_candidate_parameters() == slope_reversal_candidate_parameters()


def test_line_spectrum_returns_sorted_two_band_values():
    params = D3TwoBandParameters()
    k1 = np.linspace(-np.pi, np.pi, 23)
    k2 = np.zeros_like(k1)

    bands = line_spectrum(k1, k2, params)

    assert bands.shape == (23, 2)
    assert np.all(bands[:, 0] <= bands[:, 1])


def test_relative_filling_grid_is_centered_on_half_filling():
    grid = relative_filling_grid(-2.0, 2.0, 9)

    assert np.isclose(grid[0], -2.0)
    assert np.isclose(grid[-1], 2.0)
    assert np.isclose(grid[len(grid) // 2], 0.0)


def test_background_candidate_moves_upper_saddle_close_to_dirac_point():
    report = candidate_report(background_candidate_parameters(), nk=81)
    nearest = report["upper_saddles_above_dirac"][0]

    assert 0.03 < nearest["energy_minus_dirac"] < 0.045
    assert abs(nearest["hessian_det"]) < 0.25


def test_background_candidate_preserves_starting_line_dirac_locations():
    starting_roots = find_line_dirac_roots(starting_candidate_parameters())
    background_roots = find_line_dirac_roots(background_candidate_parameters())

    assert np.allclose(background_roots, starting_roots, atol=1e-10)


def test_candidate_report_includes_lower_saddle_isolation_diagnostics():
    report = candidate_report(background_candidate_parameters(), nk=61)

    assert "lower_saddles_below_dirac" in report
    assert report["lower_saddles_below_dirac"]
    assert report["lower_saddles_below_dirac"][0]["dirac_minus_energy"] > 0.25


def test_strong_bending_candidate_has_flat_upper_saddle_and_isolated_hole_side():
    report = candidate_report(strong_bending_candidate_parameters(), nk=61)
    nearest_upper = report["upper_saddles_above_dirac"][0]
    nearest_lower = report["lower_saddles_below_dirac"][0]

    assert 0.025 < nearest_upper["energy_minus_dirac"] < 0.045
    assert abs(nearest_upper["hessian_det"]) < 0.008
    assert nearest_lower["dirac_minus_energy"] > 0.25


def test_strong_bending_candidate_preserves_starting_line_dirac_locations():
    starting_roots = find_line_dirac_roots(starting_candidate_parameters())
    strong_roots = find_line_dirac_roots(strong_bending_candidate_parameters())

    assert np.allclose(strong_roots, starting_roots, atol=1e-10)


def test_slope_reversal_candidate_has_lower_lifshitz_scale_and_isolated_hole_side():
    report = candidate_report(slope_reversal_candidate_parameters(), nk=61)
    nearest_upper = report["upper_saddles_above_dirac"][0]
    nearest_lower = report["lower_saddles_below_dirac"][0]

    assert 0.020 < nearest_upper["energy_minus_dirac"] < 0.024
    assert abs(nearest_upper["hessian_det"]) < 0.003
    assert nearest_lower["dirac_minus_energy"] > 0.25


def test_slope_reversal_candidate_preserves_starting_line_dirac_locations():
    starting_roots = find_line_dirac_roots(starting_candidate_parameters())
    reversal_roots = find_line_dirac_roots(slope_reversal_candidate_parameters())

    assert np.allclose(reversal_roots, starting_roots, atol=1e-10)
