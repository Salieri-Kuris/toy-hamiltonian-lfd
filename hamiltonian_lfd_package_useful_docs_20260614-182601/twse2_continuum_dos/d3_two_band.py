"""D3-symmetric two-band tight-binding Hamiltonian diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .peierls_hofstadter import TBModel

S0 = np.eye(2, dtype=complex)
SX = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
SZ = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


@dataclass(frozen=True)
class D3TwoBandParameters:
    """Parameters for the periodic two-band D3 ansatz."""

    m0: float = -0.1
    m1: float = 0.2
    m2: float = -0.05
    m3: float = 0.0
    lam: float = 0.45
    u0: float = 0.0
    u1: float = -0.14072753558865014
    u2: float = -0.02245407803850541
    u3: float = -0.013573216068830357
    u4: float = 0.0335095668835709
    u5: float = 0.0
    u6: float = 0.0
    u7: float = 0.0
    internal_degeneracy: int = 2
    copy_degeneracy: int = 3


def starting_candidate_parameters() -> D3TwoBandParameters:
    """Return the current saved starting Hamiltonian candidate."""
    return D3TwoBandParameters()


def baseline_candidate_parameters() -> D3TwoBandParameters:
    """Return the earlier minimal Dirac-cone toy baseline."""
    return D3TwoBandParameters(
        u1=0.08,
        u2=-0.04,
        u3=0.0,
        u4=0.0,
        u5=0.0,
        u6=0.0,
        u7=0.0,
    )


def background_candidate_parameters() -> D3TwoBandParameters:
    """Return a tuned candidate with delayed, flatter particle-side beta DOS."""
    return D3TwoBandParameters(
        u1=-0.13786037398264536,
        u2=-0.007820340415215885,
        u3=-0.01007068900254546,
        u4=0.021857829860253575,
    )


def strong_bending_candidate_parameters() -> D3TwoBandParameters:
    """Return a high-asymmetry candidate with a very flat particle-side saddle."""
    return D3TwoBandParameters(
        u1=-0.11823260573457668,
        u2=-0.0393635379740415,
        u3=-0.023103272118251236,
        u4=-0.018348877905879905,
        u5=0.010152092446889431,
        u6=0.014070201232323309,
        u7=0.0013351112620801405,
    )


def slope_reversal_candidate_parameters() -> D3TwoBandParameters:
    """Return a more aggressive candidate tuned toward particle-side slope reversal."""
    return D3TwoBandParameters(
        u1=-0.11824074424096692,
        u2=-0.03936296247689948,
        u3=-0.023118973634458233,
        u4=-0.018372676691527402,
        u5=0.010154034407805107,
        u6=0.014014999412566741,
        u7=0.0012894556500397858,
    )


def c1_harmonic(k1: float | np.ndarray, k2: float | np.ndarray) -> np.ndarray:
    """Return cos k1 + cos k2 + cos(k1+k2)."""
    return np.cos(k1) + np.cos(k2) + np.cos(np.asarray(k1) + np.asarray(k2))


def c2_harmonic(k1: float | np.ndarray, k2: float | np.ndarray) -> np.ndarray:
    """Return the second D3 cosine harmonic."""
    return np.cos(2.0 * np.asarray(k1)) + np.cos(2.0 * np.asarray(k2)) + np.cos(
        2.0 * (np.asarray(k1) + np.asarray(k2))
    )


def g_harmonic(k1: float | np.ndarray, k2: float | np.ndarray) -> np.ndarray:
    """Return the C2'T-line-vanishing red-line sine harmonic."""
    return np.sin(np.asarray(k1) + np.asarray(k2))


def m_term(k1: float | np.ndarray, k2: float | np.ndarray, params: D3TwoBandParameters) -> np.ndarray:
    """Return M(k), the sigma_z mass term."""
    c1 = c1_harmonic(k1, k2)
    c2 = c2_harmonic(k1, k2)
    return params.m0 + params.m1 * c1 + params.m2 * c2 + params.m3 * c1 * c1


def d0_term(k1: float | np.ndarray, k2: float | np.ndarray, params: D3TwoBandParameters) -> np.ndarray:
    """Return d0(k), the scalar particle-hole-asymmetry term."""
    c1 = c1_harmonic(k1, k2)
    c2 = c2_harmonic(k1, k2)
    return (
        params.u0
        + params.u1 * c1
        + params.u2 * c2
        + params.u3 * c1 * c1
        + params.u4 * c1 * c2
        + params.u5 * c2 * c2
        + params.u6 * c1 * c1 * c1
        + params.u7 * c1 * c1 * c2
    )


def hamiltonian(k1: float, k2: float, params: D3TwoBandParameters | None = None) -> np.ndarray:
    """Return H(k)=d0*s0+M*sz+G*sx."""
    p = params or D3TwoBandParameters()
    return d0_term(k1, k2, p) * S0 + m_term(k1, k2, p) * SZ + p.lam * g_harmonic(k1, k2) * SX


def eigvals(k1: float, k2: float, params: D3TwoBandParameters | None = None) -> np.ndarray:
    """Return sorted two-band eigenvalues."""
    return np.linalg.eigvalsh(hamiltonian(k1, k2, params))


def line_spectrum(
    k1_values: np.ndarray,
    k2_values: np.ndarray,
    params: D3TwoBandParameters | None = None,
) -> np.ndarray:
    """Evaluate sorted two-band eigenvalues along a momentum path."""
    p = params or D3TwoBandParameters()
    k1_arr = np.asarray(k1_values, dtype=float)
    k2_arr = np.asarray(k2_values, dtype=float)
    if k1_arr.shape != k2_arr.shape:
        raise ValueError("k1_values and k2_values must have the same shape")
    bands = np.empty((k1_arr.size, 2), dtype=float)
    for index, (k1, k2) in enumerate(zip(k1_arr.ravel(), k2_arr.ravel())):
        bands[index] = eigvals(float(k1), float(k2), p)
    return bands


def relative_filling_grid(n_min: float, n_max: float, n_points: int) -> np.ndarray:
    """Return a relative filling grid centered on the Dirac half-filling point."""
    if n_points < 2:
        raise ValueError("n_points must be at least 2")
    if n_min >= n_max:
        raise ValueError("n_min must be smaller than n_max")
    return np.linspace(float(n_min), float(n_max), int(n_points))


FourierDict = dict[tuple[int, int], complex]


def _add_scaled(target: FourierDict, source: FourierDict, scale: complex) -> None:
    for key, value in source.items():
        target[key] = target.get(key, 0.0j) + scale * value


def _convolve(left: FourierDict, right: FourierDict) -> FourierDict:
    out: FourierDict = {}
    for (l1, l2), lvalue in left.items():
        for (r1, r2), rvalue in right.items():
            key = (l1 + r1, l2 + r2)
            out[key] = out.get(key, 0.0j) + lvalue * rvalue
    return out


def _constant() -> FourierDict:
    return {(0, 0): 1.0 + 0.0j}


def _c1_coeffs() -> FourierDict:
    coeffs: FourierDict = {}
    for key in [(1, 0), (0, 1), (1, 1)]:
        coeffs[key] = coeffs.get(key, 0.0j) + 0.5
        coeffs[(-key[0], -key[1])] = coeffs.get((-key[0], -key[1]), 0.0j) + 0.5
    return coeffs


def _c2_coeffs() -> FourierDict:
    coeffs: FourierDict = {}
    for key in [(2, 0), (0, 2), (2, 2)]:
        coeffs[key] = coeffs.get(key, 0.0j) + 0.5
        coeffs[(-key[0], -key[1])] = coeffs.get((-key[0], -key[1]), 0.0j) + 0.5
    return coeffs


def _g_coeffs() -> FourierDict:
    coeffs: FourierDict = {}
    key = (1, 1)
    coeffs[key] = coeffs.get(key, 0.0j) + 1.0 / (2.0j)
    coeffs[(-key[0], -key[1])] = coeffs.get((-key[0], -key[1]), 0.0j) - 1.0 / (2.0j)
    return coeffs


def scalar_fourier_coefficients(params: D3TwoBandParameters) -> tuple[FourierDict, FourierDict, FourierDict]:
    """Return Fourier coefficients for d0, M, and G."""
    c1 = _c1_coeffs()
    c2 = _c2_coeffs()
    c1_sq = _convolve(c1, c1)
    c1_c2 = _convolve(c1, c2)
    c2_sq = _convolve(c2, c2)
    c1_cubed = _convolve(c1_sq, c1)
    c1_sq_c2 = _convolve(c1_sq, c2)

    d0_coeffs: FourierDict = {}
    m_coeffs: FourierDict = {}
    g_coeffs: FourierDict = {}

    _add_scaled(d0_coeffs, _constant(), params.u0)
    _add_scaled(d0_coeffs, c1, params.u1)
    _add_scaled(d0_coeffs, c2, params.u2)
    _add_scaled(d0_coeffs, c1_sq, params.u3)
    _add_scaled(d0_coeffs, c1_c2, params.u4)
    _add_scaled(d0_coeffs, c2_sq, params.u5)
    _add_scaled(d0_coeffs, c1_cubed, params.u6)
    _add_scaled(d0_coeffs, c1_sq_c2, params.u7)

    _add_scaled(m_coeffs, _constant(), params.m0)
    _add_scaled(m_coeffs, c1, params.m1)
    _add_scaled(m_coeffs, c2, params.m2)
    _add_scaled(m_coeffs, c1_sq, params.m3)

    _add_scaled(g_coeffs, _g_coeffs(), params.lam)
    return d0_coeffs, m_coeffs, g_coeffs


def tb_model(params: D3TwoBandParameters | None = None, *, atol: float = 1e-14) -> TBModel:
    """Return the exact finite hopping model for the D3 two-band ansatz."""
    p = params or D3TwoBandParameters()
    d0_coeffs, m_coeffs, g_coeffs = scalar_fourier_coefficients(p)
    keys = set(d0_coeffs) | set(m_coeffs) | set(g_coeffs)
    hoppings: list[tuple[int, int, np.ndarray]] = []
    for dm, dn in sorted(keys):
        mat = d0_coeffs.get((dm, dn), 0.0j) * S0
        mat = mat + m_coeffs.get((dm, dn), 0.0j) * SZ
        mat = mat + g_coeffs.get((dm, dn), 0.0j) * SX
        if np.linalg.norm(mat) > atol:
            hoppings.append((dm, dn, np.asarray(mat, dtype=complex)))
    return TBModel(n_orb=2, hoppings=hoppings)


def _wrap_angle(value: float) -> float:
    return float((value + np.pi) % (2.0 * np.pi) - np.pi)


def _bisect_root(func, left: float, right: float, *, xtol: float = 1e-12, max_iter: int = 100) -> float:
    f_left = float(func(left))
    f_right = float(func(right))
    if abs(f_left) < xtol:
        return float(left)
    if abs(f_right) < xtol:
        return float(right)
    if f_left * f_right > 0.0:
        raise ValueError("interval does not bracket a root")
    lo = float(left)
    hi = float(right)
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = float(func(mid))
        if abs(f_mid) < xtol or abs(hi - lo) < xtol:
            return mid
        if f_left * f_mid <= 0.0:
            hi = mid
            f_right = f_mid
        else:
            lo = mid
            f_left = f_mid
    return 0.5 * (lo + hi)


def find_line_dirac_roots(
    params: D3TwoBandParameters | None = None,
    *,
    n_grid: int = 1201,
    line: str = "c2t_kprime_gamma_k",
) -> list[float]:
    """Locate roots of M(k)=0 on a C2'T-invariant red line.

    ``line="c2t_kprime_gamma_k"`` is the K-Gamma-K' line with ``k2=-k1``.
    ``line="c2t_boundary"`` is the equivalent K'-M-K boundary with
    ``k2=pi-k1``.
    """
    p = params or D3TwoBandParameters()
    grid = np.linspace(-np.pi, np.pi, n_grid)
    if line == "c2t_kprime_gamma_k":
        k2_of_k1 = lambda k1: -k1
    elif line == "c2t_boundary":
        k2_of_k1 = lambda k1: np.pi - k1
    else:
        raise ValueError(f"unknown C2'T line: {line}")
    values = np.asarray([m_term(float(k1), float(k2_of_k1(float(k1))), p) for k1 in grid], dtype=float)
    roots: list[float] = []
    func = lambda k1: m_term(k1, float(k2_of_k1(k1)), p)
    for i in range(n_grid - 1):
        left = float(grid[i])
        right = float(grid[i + 1])
        f_left = values[i]
        f_right = values[i + 1]
        if abs(f_left) < 1e-10:
            roots.append(_wrap_angle(left))
        elif f_left * f_right < 0.0:
            roots.append(_wrap_angle(_bisect_root(func, left, right)))
        elif abs(f_right) < 1e-10:
            roots.append(_wrap_angle(right))

    unique: list[float] = []
    for root in sorted(roots):
        if not unique or abs(root - unique[-1]) > 1e-5:
            unique.append(root)
    if len(unique) > 1 and abs((unique[0] + 2.0 * np.pi) - unique[-1]) < 1e-5:
        unique.pop()
    return unique


def dirac_diagnostic(
    k1: float,
    params: D3TwoBandParameters | None = None,
    *,
    line: str = "c2t_kprime_gamma_k",
    h: float = 1e-5,
) -> dict[str, float]:
    """Return energy and linear velocities at a line Dirac root."""
    p = params or D3TwoBandParameters()
    if line == "c2t_kprime_gamma_k":
        k2 = -float(k1)
        tangent = np.array([1.0, -1.0], dtype=float)
        normal = np.array([1.0, 1.0], dtype=float)
    elif line == "c2t_boundary":
        k2 = float(np.pi - k1)
        tangent = np.array([1.0, -1.0], dtype=float)
        normal = np.array([1.0, 1.0], dtype=float)
    else:
        raise ValueError(f"unknown C2'T line: {line}")
    tangent = tangent / float(np.linalg.norm(tangent))
    normal = normal / float(np.linalg.norm(normal))
    energy = float(d0_term(k1, k2, p))
    v_parallel = float(
        (
            m_term(k1 + h * tangent[0], k2 + h * tangent[1], p)
            - m_term(k1 - h * tangent[0], k2 - h * tangent[1], p)
        )
        / (2.0 * h)
    )
    v_perp = float(
        p.lam
        * (
            g_harmonic(k1 + h * normal[0], k2 + h * normal[1])
            - g_harmonic(k1 - h * normal[0], k2 - h * normal[1])
        )
        / (2.0 * h)
    )
    return {
        "k1": float(k1),
        "k2": float(_wrap_angle(k2)),
        "line": line,
        "energy": energy,
        "v_parallel": v_parallel,
        "v_perp": v_perp,
    }


def _energy_mesh(params: D3TwoBandParameters, nk: int) -> np.ndarray:
    k_values = (np.arange(nk) + 0.5) * (2.0 * np.pi / nk) - np.pi
    energies = np.empty((nk * nk, 2), dtype=float)
    index = 0
    for k1 in k_values:
        for k2 in k_values:
            energies[index] = eigvals(float(k1), float(k2), params)
            index += 1
    return energies


def filling_at_energy(
    params: D3TwoBandParameters | None,
    energy: float,
    *,
    nu_offset: float = 0.0,
    nk: int = 151,
) -> float:
    """Return active-model filling with copy/internal degeneracy weights."""
    p = params or D3TwoBandParameters()
    energies = _energy_mesh(p, nk)
    weight = p.copy_degeneracy * p.internal_degeneracy
    occupied_per_cell = weight * float(np.mean(energies <= energy) * energies.shape[1])
    return float(nu_offset + occupied_per_cell)


def calibrate_nu_offset(
    params: D3TwoBandParameters | None,
    energy: float,
    *,
    target_nu: float = -2.0,
    nk: int = 151,
) -> float:
    """Choose nu_offset so the selected energy maps to target_nu."""
    return float(target_nu - filling_at_energy(params, energy, nu_offset=0.0, nk=nk))


def _band_energy(k: np.ndarray, params: D3TwoBandParameters, band: int) -> float:
    return float(eigvals(_wrap_angle(float(k[0])), _wrap_angle(float(k[1])), params)[band])


def _gradient(k: np.ndarray, params: D3TwoBandParameters, band: int, h: float = 1e-4) -> np.ndarray:
    k1, k2 = float(k[0]), float(k[1])
    dx = (_band_energy(np.array([k1 + h, k2]), params, band) - _band_energy(np.array([k1 - h, k2]), params, band)) / (
        2.0 * h
    )
    dy = (_band_energy(np.array([k1, k2 + h]), params, band) - _band_energy(np.array([k1, k2 - h]), params, band)) / (
        2.0 * h
    )
    return np.array([dx, dy], dtype=float)


def _hessian(k: np.ndarray, params: D3TwoBandParameters, band: int, h: float = 2e-4) -> np.ndarray:
    k1, k2 = float(k[0]), float(k[1])
    f00 = _band_energy(np.array([k1, k2]), params, band)
    fpp = _band_energy(np.array([k1 + h, k2 + h]), params, band)
    fpm = _band_energy(np.array([k1 + h, k2 - h]), params, band)
    fmp = _band_energy(np.array([k1 - h, k2 + h]), params, band)
    fmm = _band_energy(np.array([k1 - h, k2 - h]), params, band)
    f_xx = (_band_energy(np.array([k1 + h, k2]), params, band) - 2.0 * f00 + _band_energy(np.array([k1 - h, k2]), params, band)) / (h * h)
    f_yy = (_band_energy(np.array([k1, k2 + h]), params, band) - 2.0 * f00 + _band_energy(np.array([k1, k2 - h]), params, band)) / (h * h)
    f_xy = (fpp - fpm - fmp + fmm) / (4.0 * h * h)
    return np.array([[f_xx, f_xy], [f_xy, f_yy]], dtype=float)


def stationary_points(
    params: D3TwoBandParameters | None = None,
    *,
    band: int = 1,
    n_seed: int = 9,
    grad_tol: float = 2e-5,
) -> list[dict[str, float]]:
    """Find coarse stationary points for a selected zero-field band."""
    p = params or D3TwoBandParameters()
    seeds = np.linspace(-np.pi, np.pi, n_seed, endpoint=False)
    points: list[dict[str, float]] = []
    for k1 in seeds:
        for k2 in seeds:
            objective = lambda x: float(np.dot(_gradient(np.asarray(x), p, band), _gradient(np.asarray(x), p, band)))
            result = minimize(objective, np.array([k1, k2], dtype=float), method="Nelder-Mead", options={"maxiter": 300})
            k = np.array([_wrap_angle(float(result.x[0])), _wrap_angle(float(result.x[1]))], dtype=float)
            grad_norm = float(np.linalg.norm(_gradient(k, p, band)))
            if grad_norm > grad_tol:
                continue
            if any(np.linalg.norm(np.angle(np.exp(1.0j * (k - np.array([q["k1"], q["k2"]]))))) < 1e-3 for q in points):
                continue
            hess = _hessian(k, p, band)
            points.append(
                {
                    "k1": float(k[0]),
                    "k2": float(k[1]),
                    "energy": _band_energy(k, p, band),
                    "grad_norm": grad_norm,
                    "hessian_det": float(np.linalg.det(hess)),
                    "hessian_trace": float(np.trace(hess)),
                }
            )
    return sorted(points, key=lambda item: item["energy"])


def candidate_report(params: D3TwoBandParameters | None = None, *, nk: int = 151) -> dict:
    """Return compact diagnostics for a D3 two-band candidate."""
    p = params or D3TwoBandParameters()
    roots = find_line_dirac_roots(p)
    if not roots:
        raise ValueError("no line Dirac roots found")
    diracs = [dirac_diagnostic(root, p, line="c2t_kprime_gamma_k") for root in roots]
    selected = min(diracs, key=lambda item: abs(item["energy"]))
    offset = calibrate_nu_offset(p, selected["energy"], target_nu=-2.0, nk=nk)
    upper_saddles = []
    for item in stationary_points(p, band=1):
        if item["hessian_det"] < 0.0 and item["energy"] > selected["energy"]:
            enriched = dict(item)
            enriched["energy_minus_dirac"] = float(item["energy"] - selected["energy"])
            upper_saddles.append(enriched)
    upper_saddles.sort(key=lambda item: item["energy_minus_dirac"])
    lower_saddles = []
    for item in stationary_points(p, band=0):
        if item["hessian_det"] < 0.0 and item["energy"] < selected["energy"]:
            enriched = dict(item)
            enriched["dirac_minus_energy"] = float(selected["energy"] - item["energy"])
            lower_saddles.append(enriched)
    lower_saddles.sort(key=lambda item: item["dirac_minus_energy"])
    return {
        "parameters": p.__dict__,
        "line_dirac_roots": diracs,
        "selected_dirac": selected,
        "nu_offset": offset,
        "selected_dirac_filling": filling_at_energy(p, selected["energy"], nu_offset=offset, nk=nk),
        "upper_saddles_above_dirac": upper_saddles,
        "lower_saddles_below_dirac": lower_saddles,
    }
