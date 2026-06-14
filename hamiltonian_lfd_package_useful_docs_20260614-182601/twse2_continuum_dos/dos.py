"""DOS broadening, filling, and density-axis utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.ndimage as ndi


@dataclass(frozen=True)
class FanResult:
    """DOS data on energy and density axes."""

    energy_grid: np.ndarray
    b_values: np.ndarray
    dos_EB: np.ndarray
    n_EB: np.ndarray
    n_grid: np.ndarray
    mu_nB: np.ndarray
    dos_nB: np.ndarray


def _gaussian(grid: np.ndarray, center: float, width: float) -> np.ndarray:
    return np.exp(-0.5 * ((grid - center) / width) ** 2) / (np.sqrt(2.0 * np.pi) * width)


def gaussian_dos(
    energy_grid: np.ndarray,
    energies: np.ndarray,
    eta: float,
    degeneracy: float = 1.0,
) -> np.ndarray:
    """Gaussian-broaden eigenvalues into DOS per moiré cell per meV."""
    grid = np.asarray(energy_grid, dtype=float)
    values = np.asarray(energies, dtype=float)
    if values.ndim == 0:
        values = values.reshape(1)

    n_k = values.shape[0] if values.ndim >= 2 else 1
    dos = np.zeros_like(grid)
    for energy in values.reshape(-1):
        dos += _gaussian(grid, float(energy), eta)
    return degeneracy * dos / n_k


def histogram_gaussian_dos(
    energy_grid: np.ndarray,
    energies: np.ndarray,
    eta: float,
    weights: np.ndarray | None = None,
    degeneracy: float = 1.0,
    truncate: float = 6.0,
) -> np.ndarray:
    """Fast weighted Gaussian DOS on a uniform energy grid.

    This bins eigenvalues to the nearest energy-grid histogram and applies a
    Gaussian filter. It is much faster than adding one continuous Gaussian per
    eigenvalue and is accurate when the grid spacing is small compared with
    ``eta``.
    """
    grid = np.asarray(energy_grid, dtype=float)
    if grid.ndim != 1 or len(grid) < 2:
        raise ValueError("energy_grid must be a one-dimensional array with at least two points")
    spacing = np.diff(grid)
    dx = float(spacing[0])
    if dx <= 0.0 or not np.allclose(spacing, dx, rtol=1e-6, atol=1e-12):
        raise ValueError("energy_grid must be uniformly spaced and increasing")
    if eta <= 0.0:
        raise ValueError("eta must be positive")

    values = np.asarray(energies, dtype=float).reshape(-1)
    finite = np.isfinite(values)
    values = values[finite]
    if weights is None:
        hist_weights = np.ones_like(values, dtype=float)
    else:
        hist_weights = np.asarray(weights, dtype=float).reshape(-1)[finite]

    edges = np.empty(len(grid) + 1, dtype=float)
    edges[1:-1] = 0.5 * (grid[:-1] + grid[1:])
    edges[0] = grid[0] - 0.5 * dx
    edges[-1] = grid[-1] + 0.5 * dx
    hist, _ = np.histogram(values, bins=edges, weights=hist_weights)
    sigma_bins = eta / dx
    broadened = ndi.gaussian_filter1d(hist.astype(float), sigma=sigma_bins, mode="constant", truncate=truncate)
    return degeneracy * broadened / dx


def saddle_window_gaussian_dos(
    energy_grid: np.ndarray,
    energies: np.ndarray,
    *,
    saddle_energy: float,
    gamma_alpha: float = 0.003,
    gamma_beta: float = 0.010,
    saddle_width: float = 0.025,
    weights: np.ndarray | None = None,
    degeneracy: float = 1.0,
    gamma_bins: int = 9,
) -> np.ndarray:
    """Gaussian DOS with larger broadening for states near a saddle energy.

    The width assigned to each eigenvalue is
    ``gamma_alpha + (gamma_beta - gamma_alpha) * exp(-0.5*((E-Es)/w)^2)``.
    Eigenvalues are grouped into a small number of width bins so this remains
    a histogram/filter calculation rather than a per-state Gaussian sum.
    """
    if gamma_alpha <= 0.0 or gamma_beta <= 0.0:
        raise ValueError("gamma_alpha and gamma_beta must be positive")
    if saddle_width <= 0.0:
        raise ValueError("saddle_width must be positive")
    if gamma_bins < 1:
        raise ValueError("gamma_bins must be at least 1")

    values = np.asarray(energies, dtype=float).reshape(-1)
    finite = np.isfinite(values)
    values = values[finite]
    if weights is None:
        state_weights = np.ones_like(values, dtype=float)
    else:
        weights_arr = np.asarray(weights, dtype=float).reshape(-1)
        if weights_arr.shape != np.asarray(energies, dtype=float).reshape(-1).shape:
            raise ValueError("weights must have the same size as energies")
        state_weights = weights_arr[finite]

    if values.size == 0:
        return np.zeros_like(np.asarray(energy_grid, dtype=float))

    widths = gamma_alpha + (gamma_beta - gamma_alpha) * np.exp(
        -0.5 * ((values - float(saddle_energy)) / saddle_width) ** 2
    )
    if gamma_bins == 1 or np.allclose(widths, widths[0], rtol=1e-12, atol=1e-15):
        return histogram_gaussian_dos(
            energy_grid,
            values,
            eta=float(np.mean(widths)),
            weights=state_weights,
            degeneracy=degeneracy,
        )

    width_min = float(np.min(widths))
    width_max = float(np.max(widths))
    edges = np.linspace(width_min, width_max, int(gamma_bins) + 1)
    dos = np.zeros_like(np.asarray(energy_grid, dtype=float))
    for index in range(int(gamma_bins)):
        if index == int(gamma_bins) - 1:
            mask = (widths >= edges[index]) & (widths <= edges[index + 1])
        else:
            mask = (widths >= edges[index]) & (widths < edges[index + 1])
        if not np.any(mask):
            continue
        eta_weights = np.abs(state_weights[mask])
        if np.any(eta_weights > 0.0):
            eta = float(np.average(widths[mask], weights=eta_weights))
        else:
            eta = float(np.mean(widths[mask]))
        dos += histogram_gaussian_dos(
            energy_grid,
            values[mask],
            eta=eta,
            weights=state_weights[mask],
            degeneracy=degeneracy,
        )
    return dos


def beta_reservoir_dos(
    energy_grid: np.ndarray,
    b_values: np.ndarray,
    *,
    center_energy: float,
    width: float = 0.08,
    amplitude: float = 0.5,
    center_slope: float = 0.0,
    center_reference_field: float = 0.0,
    b_star: float = 5.5,
    power: float = 2.0,
    onset_field: float | None = None,
    onset_width: float = 0.5,
    lower_energy: float | None = None,
    lower_gate_width: float = 0.01,
) -> np.ndarray:
    """Return a smooth particle-side beta reservoir DOS.

    The reservoir is the energy derivative of
    ``N_beta(mu, B) = amplitude * P(B) * sigmoid((mu - E_beta(B)) / width)``
    where ``E_beta(B) = center_energy + center_slope * (B - center_reference_field)``.
    with ``P(B) = 1 - exp(-(B / b_star) ** power)``. If ``onset_field`` is set,
    ``P(B)`` is replaced by a logistic switch centered at that field so the
    finite-field onset has a direct field scale. Its integral over a wide
    energy window is the field-dependent extra filling that bends dark lines.
    """
    grid = np.asarray(energy_grid, dtype=float)
    fields = np.asarray(b_values, dtype=float)
    if grid.ndim != 1 or len(grid) < 2:
        raise ValueError("energy_grid must be a one-dimensional array with at least two points")
    if width <= 0.0:
        raise ValueError("width must be positive")
    if amplitude < 0.0:
        raise ValueError("amplitude must be non-negative")
    if b_star <= 0.0:
        raise ValueError("b_star must be positive")
    if power <= 0.0:
        raise ValueError("power must be positive")
    if onset_width <= 0.0:
        raise ValueError("onset_width must be positive")
    if lower_gate_width <= 0.0:
        raise ValueError("lower_gate_width must be positive")

    centers = float(center_energy) + float(center_slope) * (fields - float(center_reference_field))
    z = (grid[None, :] - centers[:, None]) / float(width)
    sigmoid = 0.5 * (1.0 + np.tanh(0.5 * z))
    density_profile = sigmoid * (1.0 - sigmoid) / float(width)
    if lower_energy is not None:
        gate_z = (grid - float(lower_energy)) / float(lower_gate_width)
        gate = 0.5 * (1.0 + np.tanh(0.5 * gate_z))
        density_profile = density_profile * gate[None, :]
        area = np.trapz(density_profile, grid, axis=1)
        density_profile = np.divide(
            density_profile,
            area[:, None],
            out=np.zeros_like(density_profile),
            where=area[:, None] > 0.0,
        )
    if onset_field is not None:
        onset = 0.5 * (1.0 + np.tanh(0.5 * (fields - float(onset_field)) / float(onset_width)))
    else:
        onset = 1.0 - np.exp(-np.power(np.maximum(fields, 0.0) / float(b_star), float(power)))
    return float(amplitude) * onset[:, None] * density_profile


def apply_beta_hybridization_gap(
    energy_grid: np.ndarray,
    b_values: np.ndarray,
    dos_EB: np.ndarray,
    *,
    center_energy: float,
    width: float = 0.03,
    depth: float = 0.6,
    onset_field: float = 4.5,
    onset_width: float = 0.4,
    center_slope: float = 0.0,
    center_reference_field: float = 0.0,
) -> np.ndarray:
    """Suppress DOS near a finite-field beta hybridization gap.

    This is a phenomenological avoided-crossing notch in ``D(E,B)``. It is not
    a hand-drawn line in density space: the notch is defined on the energy axis
    and then passes through the same ``D(E,B) -> D(n,B)`` regridding pipeline.
    """
    grid = np.asarray(energy_grid, dtype=float)
    fields = np.asarray(b_values, dtype=float)
    dos = np.asarray(dos_EB, dtype=float)
    if grid.ndim != 1 or fields.ndim != 1 or dos.shape != (len(fields), len(grid)):
        raise ValueError("dos_EB must have shape (len(b_values), len(energy_grid))")
    if width <= 0.0:
        raise ValueError("width must be positive")
    if not (0.0 <= depth <= 1.0):
        raise ValueError("depth must be between 0 and 1")
    if onset_width <= 0.0:
        raise ValueError("onset_width must be positive")

    centers = float(center_energy) + float(center_slope) * (fields - float(center_reference_field))
    onset = 0.5 * (1.0 + np.tanh(0.5 * (fields - float(onset_field)) / float(onset_width)))
    notch = np.exp(-0.5 * ((grid[None, :] - centers[:, None]) / float(width)) ** 2)
    factor = 1.0 - float(depth) * onset[:, None] * notch
    return np.maximum(dos * factor, 0.0)


def magnetic_breakdown_replicated_spectrum(
    energies: np.ndarray,
    weights: np.ndarray,
    *,
    b_field: float,
    reference_energy: float,
    t_max: float = 0.04,
    onset_field: float = 4.5,
    onset_width: float = 0.5,
    particle_gate_width: float = 0.01,
    internal_degeneracy: float = 2.0,
    copy_degeneracy: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """Replicate C3 flavors and split only the particle side at finite B.

    The input spectrum is interpreted as one C3 flavor.  The output contains
    ``copy_degeneracy`` replicas.  Below ``reference_energy`` the replicas
    remain degenerate, giving unresolved hole-side weight
    ``internal_degeneracy * copy_degeneracy``.  Above ``reference_energy`` the
    replicas acquire finite-B offsets, making particle-side branches visible
    with spacing set by ``internal_degeneracy``.
    """
    values = np.asarray(energies, dtype=float).reshape(-1)
    base_weights = np.asarray(weights, dtype=float).reshape(-1)
    if base_weights.shape != values.shape:
        raise ValueError("weights must have the same shape as energies")
    if copy_degeneracy < 1:
        raise ValueError("copy_degeneracy must be positive")
    if internal_degeneracy <= 0.0:
        raise ValueError("internal_degeneracy must be positive")
    if t_max < 0.0:
        raise ValueError("t_max must be non-negative")
    if onset_width <= 0.0:
        raise ValueError("onset_width must be positive")
    if particle_gate_width <= 0.0:
        raise ValueError("particle_gate_width must be positive")

    if copy_degeneracy == 1:
        return values.copy(), base_weights * float(internal_degeneracy)

    offsets = np.linspace(-1.0, 1.0, int(copy_degeneracy), dtype=float)
    onset = 0.5 * (1.0 + np.tanh(0.5 * (float(b_field) - float(onset_field)) / float(onset_width)))
    particle_gate = 0.5 * (1.0 + np.tanh(0.5 * (values - float(reference_energy)) / float(particle_gate_width)))
    shifts = float(t_max) * onset * particle_gate[:, None] * offsets[None, :]
    replicated = (values[:, None] + shifts).reshape(-1)
    replicated_weights = np.repeat(base_weights * float(internal_degeneracy), int(copy_degeneracy))
    return replicated, replicated_weights


def _low_dos_minima(n_grid: np.ndarray, row: np.ndarray, n_min: float, n_max: float, percentile_cut: float) -> np.ndarray:
    n = np.asarray(n_grid, dtype=float)
    values = np.asarray(row, dtype=float)
    minima = (values[1:-1] < values[:-2]) & (values[1:-1] < values[2:])
    n_mid = n[1:-1]
    threshold = np.percentile(values, percentile_cut)
    keep = minima & (n_mid >= n_min) & (n_mid <= n_max) & (values[1:-1] <= threshold)
    return n_mid[keep]


def track_low_dos_branches(
    n_grid: np.ndarray,
    b_values: np.ndarray,
    dos_nB: np.ndarray,
    *,
    seed_n_values: list[float] | tuple[float, ...] | np.ndarray,
    n_min: float = 0.0,
    n_max: float = 1.0,
    percentile_cut: float = 45.0,
    max_step: float = 0.08,
) -> list[dict[str, np.ndarray]]:
    """Track seeded low-DOS branches through a density-field map.

    Seeds are interpreted at the largest available B. Each branch is then
    followed toward lower B by nearest local minimum, which filters out
    unconnected local minima while keeping the main dark lines.
    """
    n = np.asarray(n_grid, dtype=float)
    b = np.asarray(b_values, dtype=float)
    rows = np.asarray(dos_nB, dtype=float)
    seeds = np.asarray(seed_n_values, dtype=float)
    if n.ndim != 1 or b.ndim != 1 or rows.shape != (len(b), len(n)):
        raise ValueError("dos_nB must have shape (len(b_values), len(n_grid))")
    if len(seeds) == 0:
        return []
    if max_step <= 0.0:
        raise ValueError("max_step must be positive")

    order = np.argsort(b)[::-1]
    minima_by_index = {
        int(index): _low_dos_minima(n, rows[index], n_min, n_max, percentile_cut)
        for index in order
    }
    branches: list[dict[str, np.ndarray]] = []
    for seed in seeds:
        tracked_n = np.full(len(b), np.nan, dtype=float)
        previous = float(seed)
        for index in order:
            candidates = minima_by_index[int(index)]
            if candidates.size == 0:
                continue
            distances = np.abs(candidates - previous)
            best_index = int(np.argmin(distances))
            if distances[best_index] <= max_step:
                previous = float(candidates[best_index])
                tracked_n[index] = previous
        branches.append({"n": tracked_n, "b": b.copy(), "seed_n": np.asarray(seed)})
    return branches


def _low_dos_candidates(
    n_grid: np.ndarray,
    row: np.ndarray,
    n_min: float,
    n_max: float,
    percentile_cut: float,
    candidates_per_row: int,
) -> tuple[np.ndarray, np.ndarray]:
    n = np.asarray(n_grid, dtype=float)
    values = np.asarray(row, dtype=float)
    minima = (values[1:-1] < values[:-2]) & (values[1:-1] < values[2:])
    n_mid = n[1:-1]
    row_mid = values[1:-1]
    threshold = np.percentile(values, percentile_cut)
    keep = minima & (n_mid >= n_min) & (n_mid <= n_max) & (row_mid <= threshold)
    if not np.any(keep):
        window = (n >= n_min) & (n <= n_max)
        if not np.any(window):
            return np.array([], dtype=float), np.array([], dtype=float)
        indices = np.flatnonzero(window)
        best = indices[np.argsort(values[indices])[:candidates_per_row]]
        return n[best], values[best]
    candidates_n = n_mid[keep]
    candidates_cost = row_mid[keep]
    order = np.argsort(candidates_cost)[:candidates_per_row]
    return candidates_n[order], candidates_cost[order]


def track_low_dos_branches_dp(
    n_grid: np.ndarray,
    b_values: np.ndarray,
    dos_nB: np.ndarray,
    *,
    seed_n_values: list[float] | tuple[float, ...] | np.ndarray,
    seed_slopes: list[float] | tuple[float, ...] | np.ndarray | None = None,
    n_min: float = 0.0,
    n_max: float = 1.0,
    percentile_cut: float = 55.0,
    candidates_per_row: int = 16,
    seed_weight: float = 200.0,
    continuity_weight: float = 50.0,
    slope_weight: float = 0.0,
    max_step: float | None = None,
) -> list[dict[str, np.ndarray]]:
    """Track low-DOS branches by global dynamic programming from low-field seeds.

    Seeds are interpreted at the lowest available field. The optimizer chooses
    one local DOS minimum per field row with penalties for local DOS, distance
    from the seed at low field, branch continuity, and optional low-field slope.
    """
    n = np.asarray(n_grid, dtype=float)
    b = np.asarray(b_values, dtype=float)
    rows = np.asarray(dos_nB, dtype=float)
    seeds = np.asarray(seed_n_values, dtype=float)
    if n.ndim != 1 or b.ndim != 1 or rows.shape != (len(b), len(n)):
        raise ValueError("dos_nB must have shape (len(b_values), len(n_grid))")
    if len(seeds) == 0:
        return []
    if candidates_per_row < 1:
        raise ValueError("candidates_per_row must be positive")
    if seed_weight < 0.0 or continuity_weight < 0.0 or slope_weight < 0.0:
        raise ValueError("weights must be non-negative")
    if max_step is not None and max_step <= 0.0:
        raise ValueError("max_step must be positive when set")

    if seed_slopes is None:
        slopes = np.zeros_like(seeds)
        use_slope = np.zeros_like(seeds, dtype=bool)
    else:
        slopes = np.asarray(seed_slopes, dtype=float)
        if slopes.shape != seeds.shape:
            raise ValueError("seed_slopes must have the same shape as seed_n_values")
        use_slope = np.ones_like(seeds, dtype=bool)

    order = np.argsort(b)
    b_sorted = b[order]
    row_candidates: list[tuple[np.ndarray, np.ndarray]] = []
    for index in order:
        candidates, costs = _low_dos_candidates(
            n,
            rows[int(index)],
            n_min=n_min,
            n_max=n_max,
            percentile_cut=percentile_cut,
            candidates_per_row=candidates_per_row,
        )
        if candidates.size == 0:
            candidates = np.array([np.clip(0.5 * (n_min + n_max), n_min, n_max)], dtype=float)
            costs = np.array([0.0], dtype=float)
        spread = np.percentile(rows[int(index)], 90.0) - np.percentile(rows[int(index)], 10.0)
        if spread <= 0.0:
            local_cost = np.zeros_like(costs)
        else:
            local_cost = (costs - np.min(rows[int(index)])) / spread
        row_candidates.append((candidates, local_cost))

    branches: list[dict[str, np.ndarray]] = []
    b0 = float(b_sorted[0])
    for seed, slope, has_slope in zip(seeds, slopes, use_slope):
        costs = row_candidates[0][1] + seed_weight * (row_candidates[0][0] - float(seed)) ** 2
        parents: list[np.ndarray] = []
        history_costs: list[np.ndarray] = [costs]
        for row_index in range(1, len(row_candidates)):
            previous_n, previous_cost = row_candidates[row_index - 1]
            current_n, local_cost = row_candidates[row_index]
            db = float(b_sorted[row_index] - b_sorted[row_index - 1])
            transition = continuity_weight * (current_n[:, None] - previous_n[None, :]) ** 2
            if max_step is not None:
                transition = np.where(
                    np.abs(current_n[:, None] - previous_n[None, :]) <= float(max_step),
                    transition,
                    np.inf,
                )
            if has_slope and slope_weight > 0.0 and db != 0.0:
                transition += slope_weight * ((current_n[:, None] - previous_n[None, :]) / db - float(slope)) ** 2
            total = local_cost[:, None] + history_costs[-1][None, :] + transition
            parent = np.argmin(total, axis=1)
            best = total[np.arange(len(current_n)), parent]
            if not np.any(np.isfinite(best)):
                transition = continuity_weight * (current_n[:, None] - previous_n[None, :]) ** 2
                total = local_cost[:, None] + history_costs[-1][None, :] + transition
                parent = np.argmin(total, axis=1)
                best = total[np.arange(len(current_n)), parent]
            low_field_prior = seed_weight * 0.01 * (current_n - (float(seed) + float(slope) * (b_sorted[row_index] - b0))) ** 2
            history_costs.append(best + low_field_prior)
            parents.append(parent)

        selected_indices = [int(np.argmin(history_costs[-1]))]
        for parent in reversed(parents):
            selected_indices.append(int(parent[selected_indices[-1]]))
        selected_indices.reverse()

        tracked_sorted = np.empty(len(b_sorted), dtype=float)
        for row_index, candidate_index in enumerate(selected_indices):
            tracked_sorted[row_index] = row_candidates[row_index][0][candidate_index]
        tracked = np.full(len(b), np.nan, dtype=float)
        tracked[order] = tracked_sorted
        branches.append({"n": tracked, "b": b.copy(), "seed_n": np.asarray(seed), "seed_slope": np.asarray(slope)})
    return branches


def cumulative_trapezoid(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Cumulative integral with the first value fixed to zero."""
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    increments = 0.5 * (y_arr[1:] + y_arr[:-1]) * np.diff(x_arr)
    return np.concatenate(([0.0], np.cumsum(increments)))


def regrid_dos_to_density(
    energy_grid: np.ndarray,
    b_values: np.ndarray,
    dos_EB: np.ndarray,
    n_grid: np.ndarray,
    reference_energy: float | None = None,
) -> FanResult:
    """Convert D(E,B) to D(n,B) by inverting n(E,B) independently at each B."""
    energy = np.asarray(energy_grid, dtype=float)
    b_arr = np.asarray(b_values, dtype=float)
    dos = np.asarray(dos_EB, dtype=float)
    density = np.asarray(n_grid, dtype=float)

    n_rows = []
    for row in dos:
        cumulative = cumulative_trapezoid(energy, row)
        if reference_energy is not None:
            cumulative = cumulative - np.interp(reference_energy, energy, cumulative)
        n_rows.append(cumulative)
    n_EB = np.vstack(n_rows)
    mu_nB = np.zeros((len(b_arr), len(density)), dtype=float)
    dos_nB = np.zeros_like(mu_nB)

    for i, row in enumerate(dos):
        mu_nB[i] = np.interp(density, n_EB[i], energy)
        dos_nB[i] = np.interp(mu_nB[i], energy, row)

    return FanResult(
        energy_grid=energy,
        b_values=b_arr,
        dos_EB=dos,
        n_EB=n_EB,
        n_grid=density,
        mu_nB=mu_nB,
        dos_nB=dos_nB,
    )
