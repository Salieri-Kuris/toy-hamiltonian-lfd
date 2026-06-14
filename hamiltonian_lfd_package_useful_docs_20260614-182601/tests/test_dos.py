import numpy as np

from twse2_continuum_dos.dos import (
    apply_beta_hybridization_gap,
    beta_reservoir_dos,
    cumulative_trapezoid,
    gaussian_dos,
    histogram_gaussian_dos,
    magnetic_breakdown_replicated_spectrum,
    saddle_window_gaussian_dos,
    track_low_dos_branches_dp,
    track_low_dos_branches,
    regrid_dos_to_density,
)


def test_gaussian_dos_integrates_to_number_of_states_times_degeneracy():
    energies = np.array([[0.0], [1.0], [2.0]])
    energy_grid = np.linspace(-3.0, 5.0, 2001)

    dos = gaussian_dos(energy_grid, energies, eta=0.05, degeneracy=2.0)
    integral = np.trapz(dos, energy_grid)

    assert np.isclose(integral, 2.0, atol=1e-4)


def test_histogram_gaussian_dos_integrates_weighted_states():
    energies = np.array([0.0, 1.0, 2.0])
    weights = np.array([1.0, 0.5, 0.25])
    energy_grid = np.linspace(-3.0, 5.0, 2001)

    dos = histogram_gaussian_dos(energy_grid, energies, eta=0.05, weights=weights, degeneracy=2.0)
    integral = np.trapz(dos, energy_grid)

    assert np.isclose(integral, 2.0 * weights.sum(), rtol=2e-3)


def test_histogram_gaussian_dos_matches_direct_gaussian_for_dense_grid():
    energies = np.array([-0.2, 0.15, 0.7])
    energy_grid = np.linspace(-1.5, 1.5, 3001)

    direct = gaussian_dos(energy_grid, energies, eta=0.06, degeneracy=1.0)
    hist = histogram_gaussian_dos(energy_grid, energies, eta=0.06, degeneracy=1.0)

    assert np.allclose(hist, direct, atol=0.035)


def test_saddle_window_gaussian_dos_preserves_weighted_state_count():
    energies = np.array([-0.4, -0.01, 0.0, 0.012, 0.45])
    weights = np.array([1.0, 0.5, 0.75, 0.25, 1.25])
    energy_grid = np.linspace(-1.0, 1.0, 4001)

    dos = saddle_window_gaussian_dos(
        energy_grid,
        energies,
        weights=weights,
        saddle_energy=0.0,
        gamma_alpha=0.02,
        gamma_beta=0.08,
        saddle_width=0.04,
        degeneracy=2.0,
    )
    integral = np.trapz(dos, energy_grid)

    assert np.isclose(integral, 2.0 * weights.sum(), rtol=2e-3)


def test_saddle_window_gaussian_dos_broadens_saddle_states_more_than_alpha_states():
    energy_grid = np.linspace(-0.3, 0.3, 3001)

    alpha_dos = saddle_window_gaussian_dos(
        energy_grid,
        np.array([0.20]),
        saddle_energy=0.0,
        gamma_alpha=0.01,
        gamma_beta=0.08,
        saddle_width=0.02,
    )
    beta_dos = saddle_window_gaussian_dos(
        energy_grid,
        np.array([0.0]),
        saddle_energy=0.0,
        gamma_alpha=0.01,
        gamma_beta=0.08,
        saddle_width=0.02,
    )

    assert beta_dos.max() < alpha_dos.max()


def test_beta_reservoir_dos_integrates_to_field_dependent_amplitude():
    energy_grid = np.linspace(-1.0, 1.0, 4001)
    b_values = np.array([0.0, 3.0, 6.0])

    rows = beta_reservoir_dos(
        energy_grid,
        b_values,
        center_energy=0.1,
        width=0.04,
        amplitude=0.5,
        b_star=6.0,
        power=2.0,
    )
    integrals = np.trapz(rows, energy_grid, axis=1)
    expected = 0.5 * (1.0 - np.exp(-((b_values / 6.0) ** 2.0)))

    assert rows.shape == (3, len(energy_grid))
    assert np.allclose(integrals, expected, atol=2e-4)


def test_beta_reservoir_dos_requires_positive_shape_parameters():
    energy_grid = np.linspace(-1.0, 1.0, 101)

    with np.testing.assert_raises(ValueError):
        beta_reservoir_dos(energy_grid, np.array([1.0]), center_energy=0.0, width=0.0)

    with np.testing.assert_raises(ValueError):
        beta_reservoir_dos(energy_grid, np.array([1.0]), center_energy=0.0, b_star=0.0)


def test_beta_reservoir_dos_can_be_gated_to_particle_side():
    energy_grid = np.linspace(-1.0, 1.0, 4001)
    b_values = np.array([6.0])

    rows = beta_reservoir_dos(
        energy_grid,
        b_values,
        center_energy=0.08,
        width=0.12,
        amplitude=0.7,
        b_star=5.0,
        lower_energy=0.0,
        lower_gate_width=0.003,
    )
    below_dirac_weight = np.trapz(rows[0, energy_grid < 0.0], energy_grid[energy_grid < 0.0])
    total_weight = np.trapz(rows[0], energy_grid)

    assert below_dirac_weight < 0.01 * total_weight
    assert np.isclose(total_weight, 0.7 * (1.0 - np.exp(-((6.0 / 5.0) ** 2.0))), rtol=2e-3)


def test_beta_reservoir_dos_can_have_finite_field_onset():
    energy_grid = np.linspace(-1.0, 1.0, 4001)
    b_values = np.array([0.0, 3.0, 4.5, 6.0, 9.0])

    rows = beta_reservoir_dos(
        energy_grid,
        b_values,
        center_energy=0.08,
        width=0.10,
        amplitude=0.9,
        b_star=4.5,
        onset_field=4.5,
        onset_width=0.35,
    )
    integrals = np.trapz(rows, energy_grid, axis=1)

    assert integrals[1] < 0.03 * integrals[-1]
    assert 0.35 * integrals[-1] < integrals[2] < 0.65 * integrals[-1]
    assert integrals[3] > 0.80 * integrals[-1]


def test_beta_reservoir_dos_can_shift_center_with_field():
    energy_grid = np.linspace(-0.5, 0.8, 2601)
    b_values = np.array([2.0, 6.0])

    rows = beta_reservoir_dos(
        energy_grid,
        b_values,
        center_energy=0.05,
        width=0.05,
        amplitude=1.0,
        onset_field=0.0,
        center_slope=0.025,
        center_reference_field=2.0,
    )
    peak_energies = energy_grid[np.argmax(rows, axis=1)]

    assert np.isclose(peak_energies[0], 0.05, atol=0.003)
    assert np.isclose(peak_energies[1], 0.15, atol=0.003)


def test_apply_beta_hybridization_gap_suppresses_dos_after_onset():
    energy_grid = np.linspace(-0.4, 0.6, 1001)
    b_values = np.array([2.0, 6.0])
    dos = np.ones((2, len(energy_grid)))

    suppressed = apply_beta_hybridization_gap(
        energy_grid,
        b_values,
        dos,
        center_energy=0.1,
        width=0.03,
        depth=0.7,
        onset_field=4.0,
        onset_width=0.2,
    )
    center_index = int(np.argmin(np.abs(energy_grid - 0.1)))

    assert np.isclose(suppressed[0, center_index], 1.0, atol=2e-2)
    assert suppressed[1, center_index] < 0.35
    assert np.all(suppressed >= 0.0)


def test_apply_beta_hybridization_gap_center_can_move_with_field():
    energy_grid = np.linspace(-0.4, 0.6, 1001)
    b_values = np.array([4.0, 8.0])
    dos = np.ones((2, len(energy_grid)))

    suppressed = apply_beta_hybridization_gap(
        energy_grid,
        b_values,
        dos,
        center_energy=0.0,
        width=0.025,
        depth=0.8,
        onset_field=0.0,
        center_slope=0.04,
        center_reference_field=4.0,
    )
    minima = energy_grid[np.argmin(suppressed, axis=1)]

    assert np.isclose(minima[0], 0.0, atol=0.003)
    assert np.isclose(minima[1], 0.16, atol=0.003)


def test_magnetic_breakdown_replicas_keep_hole_copies_degenerate():
    energies = np.array([-0.2])
    weights = np.array([0.5])

    replicated, replicated_weights = magnetic_breakdown_replicated_spectrum(
        energies,
        weights,
        b_field=8.0,
        reference_energy=0.0,
        t_max=0.06,
        onset_field=4.0,
        onset_width=0.3,
        particle_gate_width=0.01,
        internal_degeneracy=2.0,
        copy_degeneracy=3,
    )

    assert replicated.shape == (3,)
    assert np.allclose(replicated, -0.2, atol=1e-6)
    assert np.allclose(replicated_weights, np.full(3, 1.0))
    assert np.isclose(replicated_weights.sum(), 3.0)


def test_magnetic_breakdown_replicas_split_particle_copies_after_onset():
    energies = np.array([0.2])
    weights = np.array([0.5])

    low_field, _ = magnetic_breakdown_replicated_spectrum(
        energies,
        weights,
        b_field=1.0,
        reference_energy=0.0,
        t_max=0.06,
        onset_field=4.0,
        onset_width=0.3,
    )
    high_field, high_weights = magnetic_breakdown_replicated_spectrum(
        energies,
        weights,
        b_field=8.0,
        reference_energy=0.0,
        t_max=0.06,
        onset_field=4.0,
        onset_width=0.3,
    )

    assert np.allclose(low_field, 0.2, atol=1e-5)
    assert np.allclose(np.sort(high_field), np.array([0.14, 0.2, 0.26]), atol=2e-4)
    assert np.allclose(high_weights, np.full(3, 1.0))


def test_track_low_dos_branches_keeps_seeded_main_dark_lines():
    n_grid = np.linspace(0.0, 1.0, 501)
    b_values = np.linspace(1.0, 10.0, 37)
    branch_slopes = [0.42, 0.68]
    seed_values = [slope + 0.08 for slope in branch_slopes]
    dos_rows = []
    for b in b_values:
        t = b / b_values[-1]
        centers = [slope * t + 0.08 * t * t for slope in branch_slopes]
        row = np.ones_like(n_grid)
        for center in centers:
            row -= 0.8 * np.exp(-0.5 * ((n_grid - center) / 0.008) ** 2)
        row -= 0.15 * np.exp(-0.5 * ((n_grid - 0.93) / 0.006) ** 2)
        dos_rows.append(row)

    branches = track_low_dos_branches(
        n_grid,
        b_values,
        np.vstack(dos_rows),
        seed_n_values=seed_values,
        n_min=0.02,
        n_max=0.98,
        max_step=0.08,
    )

    assert len(branches) == 2
    for branch, seed in zip(branches, seed_values):
        finite = np.isfinite(branch["n"])
        assert np.count_nonzero(finite) > 20
        assert np.isclose(branch["n"][finite][-1], seed, atol=0.02)


def test_track_low_dos_branches_dp_tracks_low_field_seeded_curved_valley():
    n_grid = np.linspace(0.0, 1.0, 801)
    b_values = np.linspace(0.8, 10.0, 47)
    true_branch = 0.07 + 0.035 * b_values - 0.0022 * b_values**2
    dos_rows = []
    for i, center in enumerate(true_branch):
        row = np.ones_like(n_grid)
        row -= 0.8 * np.exp(-0.5 * ((n_grid - center) / 0.007) ** 2)
        jump_distractor = 0.74 - 0.018 * b_values[i]
        row -= 0.76 * np.exp(-0.5 * ((n_grid - jump_distractor) / 0.006) ** 2)
        row += 0.02 * np.sin(7.0 * n_grid + 0.3 * i)
        dos_rows.append(row)

    branches = track_low_dos_branches_dp(
        n_grid,
        b_values,
        np.vstack(dos_rows),
        seed_n_values=[true_branch[0]],
        seed_slopes=[0.035],
        n_min=0.02,
        n_max=0.58,
        percentile_cut=60.0,
        candidates_per_row=8,
        seed_weight=400.0,
        continuity_weight=80.0,
        slope_weight=20.0,
        max_step=0.06,
    )

    tracked = branches[0]["n"]

    assert np.count_nonzero(np.isfinite(tracked)) == len(b_values)
    assert np.max(np.abs(tracked - true_branch)) < 0.018


def test_cumulative_trapezoid_starts_at_zero_and_integrates_constant():
    x = np.linspace(0.0, 2.0, 11)
    y = np.ones_like(x) * 3.0

    cumulative = cumulative_trapezoid(x, y)

    assert cumulative[0] == 0.0
    assert np.isclose(cumulative[-1], 6.0)


def test_single_ladder_density_fan_has_linear_degeneracy_steps():
    energy_grid = np.linspace(-1.0, 3.0, 2001)
    b_values = np.array([1.0, 2.0])
    dos_eb = np.zeros((len(b_values), len(energy_grid)))

    for i, weight in enumerate([0.1, 0.2]):
        levels = np.array([0.0, 1.0, 2.0])
        for level in levels:
            dos_eb[i] += gaussian_dos(energy_grid, np.array([level]), eta=0.02, degeneracy=weight)

    n_grid = np.array([0.05, 0.15, 0.25])
    result = regrid_dos_to_density(energy_grid, b_values, dos_eb, n_grid)

    assert result.dos_nB.shape == (2, 3)
    assert np.all(np.isfinite(result.mu_nB))
    assert result.n_EB[1, -1] > result.n_EB[0, -1]


def test_regrid_dos_to_density_can_use_zero_energy_as_signed_density_reference():
    energy_grid = np.linspace(-2.0, 2.0, 401)
    b_values = np.array([1.0])
    dos_eb = np.ones((1, len(energy_grid)))
    n_grid = np.array([-1.0, 0.0, 1.0])

    result = regrid_dos_to_density(
        energy_grid,
        b_values,
        dos_eb,
        n_grid,
        reference_energy=0.0,
    )

    zero_index = np.argmin(np.abs(energy_grid))

    assert np.isclose(result.n_EB[0, zero_index], 0.0)
    assert np.allclose(result.mu_nB[0], n_grid)
    assert np.allclose(result.dos_nB[0], np.ones_like(n_grid))
