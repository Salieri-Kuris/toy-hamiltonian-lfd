import numpy as np

from twse2_continuum_dos.peierls_hofstadter import (
    TBModel,
    magnetic_hamiltonian,
    peierls_phase,
    rational_flux_grid_for_b_field,
    spectrum_at_flux,
)


def square_lattice_model(t: float = 1.0) -> TBModel:
    onsite = np.array([[0.0]], dtype=complex)
    tx = np.array([[t]], dtype=complex)
    ty = np.array([[t]], dtype=complex)
    return TBModel(
        n_orb=1,
        hoppings=[
            (0, 0, onsite),
            (1, 0, tx),
            (-1, 0, tx.conj().T),
            (0, 1, ty),
            (0, -1, ty.conj().T),
        ],
    )


def test_tbmodel_reconstructs_zero_field_bloch_hamiltonian():
    model = square_lattice_model(t=0.7)
    k1, k2 = 0.31, -1.2

    value = model.bloch_hamiltonian(k1, k2)[0, 0]

    assert np.allclose(value, 2.0 * 0.7 * (np.cos(k1) + np.cos(k2)))
    assert model.hermiticity_error() < 1e-12


def test_peierls_phase_has_expected_primitive_plaquette_curl():
    p, q = 2, 17
    phi = p / q
    curl = (
        peierls_phase(0, 0, 1, 0, phi)
        + peierls_phase(1, 0, 0, 1, phi)
        - peierls_phase(0, 1, 1, 0, phi)
        - peierls_phase(0, 0, 0, 1, phi)
    )

    assert np.isclose(curl, 2.0 * np.pi * phi)


def test_magnetic_hamiltonian_is_hermitian_and_has_q_orbital_dimension():
    model = square_lattice_model()

    h_mag = magnetic_hamiltonian(model, p=1, q=11, k1_mag=0.02, k2_mag=0.31, sparse=False)

    assert h_mag.shape == (11, 11)
    assert np.allclose(h_mag, h_mag.conj().T, atol=1e-10)


def test_spectrum_at_flux_returns_sorted_dense_values_and_normalized_weight():
    model = square_lattice_model()
    result = spectrum_at_flux(model, p=1, q=7, kmesh=(1, 3), method="dense_full")

    assert result.eigenvalues.shape == (3, 7)
    assert result.state_weight_per_original_cell == 1.0 / (7 * 3)
    assert np.all(np.diff(result.eigenvalues, axis=1) >= -1e-12)


def test_rational_flux_grid_approximates_uniform_b_targets():
    p_values, q_values, b_values = rational_flux_grid_for_b_field(
        b_min=0.76,
        b_max=10.0,
        n_points=80,
        q_max=263,
        flux_quantum_field_tesla=200.0,
    )

    assert p_values.shape == q_values.shape == b_values.shape
    assert len(b_values) > 60
    assert b_values[0] >= 0.75
    assert b_values[-1] <= 10.05
    assert np.all(np.diff(b_values) > 0.0)
    assert np.max(q_values) <= 263
    assert np.all([np.gcd(int(p), int(q)) == 1 for p, q in zip(p_values, q_values)])

    old_b = 200.0 / np.arange(20, 264)
    old_b = np.sort(old_b[(old_b >= b_values[0]) & (old_b <= b_values[-1])])
    assert np.std(np.diff(b_values)) < np.std(np.diff(old_b))


def test_tbmodel_rejects_non_square_hopping_matrices():
    bad = np.zeros((2, 3), dtype=complex)

    try:
        TBModel(n_orb=2, hoppings=[(0, 0, bad)])
    except ValueError as exc:
        assert "shape" in str(exc)
    else:
        raise AssertionError("TBModel accepted a non-square hopping matrix")
