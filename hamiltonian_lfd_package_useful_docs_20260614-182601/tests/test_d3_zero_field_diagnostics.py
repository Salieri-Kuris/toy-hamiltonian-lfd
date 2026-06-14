import numpy as np

from twse2_continuum_dos.d3_two_band import strong_bending_candidate_parameters
from twse2_continuum_dos.d3_zero_field_diagnostics import (
    contour_component_diagnostics,
    high_symmetry_spectrum,
    upper_band_mesh,
)


def test_high_symmetry_spectrum_returns_named_paths_and_two_bands():
    result = high_symmetry_spectrum(strong_bending_candidate_parameters(), points_per_segment=24)

    assert len(result.paths) == 2
    for path in result.paths:
        assert path.x.ndim == 1
        assert path.energies.shape == (len(path.x), 2)
        assert np.all(path.energies[:, 0] <= path.energies[:, 1])


def test_upper_band_mesh_returns_square_grid_and_upper_band():
    mesh = upper_band_mesh(strong_bending_candidate_parameters(), nk=31)

    assert mesh.k1.shape == mesh.k2.shape == mesh.upper.shape == (31, 31)
    assert np.isfinite(mesh.upper).all()


def test_contour_component_diagnostics_reports_particle_side_offsets():
    result = contour_component_diagnostics(
        strong_bending_candidate_parameters(),
        offsets=np.array([0.02, 0.03]),
        nk=61,
        shell_width=0.006,
        min_pixels=5,
    )

    assert result.energy_dirac < result.levels[0]
    assert len(result.levels) == 2
    assert len(result.by_level) == 2
    assert all(item.level > result.energy_dirac for item in result.by_level)
    assert all(item.component_count >= item.alpha_candidate_count for item in result.by_level)
