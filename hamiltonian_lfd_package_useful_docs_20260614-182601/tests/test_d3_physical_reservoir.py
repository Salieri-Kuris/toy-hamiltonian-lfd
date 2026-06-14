import numpy as np

from twse2_continuum_dos.d3_physical_reservoir import (
    PhysicalReservoirParams,
    bands,
    construct_single_reservoir_candidate,
    d0_value,
    mass_value,
    offdiag_value,
    particle_breakdown_metrics,
    reservoir_window_metrics,
    tb_model,
)


def test_constructed_model_places_dirac_on_c2t_line():
    params = construct_single_reservoir_candidate()

    assert np.isclose(params.root, 1.00)
    assert np.isclose(d0_value(params.root, -params.root, params), 0.0, atol=1e-12)
    lower, upper = bands(params.root, -params.root, params)
    assert np.isclose(lower, 0.0, atol=1e-12)
    assert np.isclose(upper, 0.0, atol=1e-12)


def test_constructed_model_has_clear_untilted_half_filling_dirac_cone():
    params = construct_single_reservoir_candidate()
    h = 1e-5

    tangent = np.array([1.0, -1.0]) / np.sqrt(2.0)
    normal = np.array([1.0, 1.0]) / np.sqrt(2.0)
    vm = (
        mass_value(params.root + h * tangent[0], -params.root + h * tangent[1], params)
        - mass_value(params.root - h * tangent[0], -params.root - h * tangent[1], params)
    ) / (2.0 * h)
    vg = (
        offdiag_value(params.root + h * normal[0], -params.root + h * normal[1], params)
        - offdiag_value(params.root - h * normal[0], -params.root - h * normal[1], params)
    ) / (2.0 * h)

    assert abs(vm) > 0.03
    assert abs(vg) > 0.03


def test_constructed_model_has_nonflat_boundary_dispersion():
    params = construct_single_reservoir_candidate()
    k1 = np.linspace(np.pi / 3.0, 2.0 * np.pi / 3.0, 121)
    k2 = np.pi - k1

    lower, upper = bands(k1, k2, params)

    assert np.ptp(lower) > 0.012
    assert np.ptp(upper) > 0.012


def test_constructed_model_places_dirac_energy_at_two_band_half_filling():
    params = construct_single_reservoir_candidate()
    grid = np.linspace(-np.pi, np.pi, 81, endpoint=False)
    k1, k2 = np.meshgrid(grid, grid, indexing="ij")

    lower, upper = bands(k1, k2, params)

    assert np.max(lower) <= 1e-10
    assert np.min(upper) >= -1e-10
    occupied = np.mean(np.concatenate([(lower <= 0.0).ravel(), (upper <= 0.0).ravel()]))
    assert np.isclose(occupied, 0.5, atol=1.0e-3)


def test_constructed_model_has_one_dominant_particle_reservoir_window():
    params = construct_single_reservoir_candidate()

    metrics = reservoir_window_metrics(params, nk=181)

    assert metrics.center_window_area > 0.005
    assert metrics.side_window_area < 0.0015
    assert metrics.max_side_shell_area < 0.0008
    assert metrics.center_to_side_ratio > 20.0


def test_constructed_model_keeps_hole_side_clean_near_dirac():
    params = construct_single_reservoir_candidate()
    breakdown = particle_breakdown_metrics(params, nk=201)

    assert breakdown.hole_background_area < 0.025
    assert breakdown.hole_max_shell_area < 0.010


def test_constructed_model_has_particle_only_breakdown_corridor():
    params = construct_single_reservoir_candidate()
    breakdown = particle_breakdown_metrics(params, nk=201)

    assert breakdown.low_particle_alpha_components == 2
    assert breakdown.low_particle_background_components == 0
    assert breakdown.reservoir_components == 1
    assert 0.020 < breakdown.neck_energy < 0.075
    assert breakdown.particle_hole_asymmetry > 2.5


def test_side_pocket_metric_detects_off_center_particle_pockets():
    params = PhysicalReservoirParams(
        root=1.2,
        reservoir_energy=0.04,
        alpha_velocity=0.03,
        offdiag_velocity=0.08,
        boundary_dispersion=0.0,
    )

    metrics = reservoir_window_metrics(params, nk=151)

    assert metrics.side_window_area > 0.002


def test_physical_reservoir_public_exports_are_available():
    import twse2_continuum_dos as tw

    params = tw.d3_construct_single_reservoir_candidate(root=1.72)
    metrics = tw.d3_reservoir_window_metrics(params, nk=101)
    breakdown = tw.d3_particle_breakdown_metrics(params, nk=101)

    assert isinstance(params, tw.PhysicalReservoirParams)
    assert isinstance(metrics, tw.ReservoirWindowMetrics)
    assert isinstance(breakdown, tw.ParticleBreakdownMetrics)


def test_tb_model_reconstructs_physical_reservoir_hamiltonian():
    params = construct_single_reservoir_candidate()
    model = tb_model(params)

    assert model.n_orb == 2
    for k1, k2 in [(0.0, 0.0), (params.root, -params.root), (0.37, -1.2)]:
        assert np.allclose(model.bloch_hamiltonian(k1, k2), tw_matrix(k1, k2, params), atol=1e-12)


def tw_matrix(k1, k2, params):
    from twse2_continuum_dos.d3_physical_reservoir import hamiltonian

    return hamiltonian(k1, k2, params)
