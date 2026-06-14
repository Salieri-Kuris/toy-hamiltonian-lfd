"""Generic Peierls-Hofstadter builder for finite-range tight-binding models."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import floor
from typing import Iterable

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass(frozen=True)
class TBModel:
    """Finite-range tight-binding model on integer lattice coordinates.

    The Bloch Hamiltonian convention is

    ``H(k) = sum_(dm,dn) h_(dm,dn) exp[i (k1 dm + k2 dn)]``.

    Each hopping matrix maps orbitals in cell ``(m,n)`` to orbitals in
    ``(m+dm,n+dn)``.
    """

    n_orb: int
    hoppings: list[tuple[int, int, np.ndarray]]
    orbital_pos: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.n_orb <= 0:
            raise ValueError("n_orb must be positive")
        for _, _, mat in self.hoppings:
            if mat.shape != (self.n_orb, self.n_orb):
                raise ValueError(f"hopping matrix has shape {mat.shape}, expected {(self.n_orb, self.n_orb)}")
        if self.orbital_pos is not None and np.asarray(self.orbital_pos).shape != (self.n_orb, 2):
            raise ValueError("orbital_pos must have shape (n_orb, 2)")

    @property
    def positions(self) -> np.ndarray:
        if self.orbital_pos is None:
            return np.zeros((self.n_orb, 2), dtype=float)
        return np.asarray(self.orbital_pos, dtype=float)

    def bloch_hamiltonian(self, k1: float, k2: float) -> np.ndarray:
        """Return the zero-field Bloch Hamiltonian."""
        out = np.zeros((self.n_orb, self.n_orb), dtype=complex)
        for dm, dn, mat in self.hoppings:
            out += mat * np.exp(1.0j * (k1 * dm + k2 * dn))
        return out

    def hermiticity_error(self, samples: Iterable[tuple[float, float]] | None = None) -> float:
        """Maximum Bloch-Hamiltonian Hermiticity residual over sample momenta."""
        if samples is None:
            samples = [(0.0, 0.0), (0.31, -0.27), (1.7, 2.1), (-2.4, 0.8)]
        errors = []
        for k1, k2 in samples:
            h = self.bloch_hamiltonian(k1, k2)
            errors.append(float(la.norm(h - h.conj().T)))
        return max(errors)


@dataclass(frozen=True)
class HofstadterSpectrum:
    """Eigenvalues and optional eigenvectors from a magnetic flux calculation."""

    p: int
    q: int
    k_points: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray | None
    state_weight_per_original_cell: float


def peierls_phase(
    m: int,
    n: int,
    dm: int,
    dn: int,
    phi: float,
    orbital_from: tuple[float, float] = (0.0, 0.0),
    orbital_to: tuple[float, float] = (0.0, 0.0),
) -> float:
    """Lattice-coordinate Landau-gauge Peierls phase for one hopping.

    The hopping goes from orbital ``from`` in cell ``(m,n)`` to orbital ``to``
    in cell ``(m+dm,n+dn)``. With all orbitals at the cell origin this reduces
    to ``2*pi*phi*(m + dm/2)*dn``.
    """
    x_beta, y_beta = orbital_from
    x_alpha, y_alpha = orbital_to
    return float(
        2.0
        * np.pi
        * phi
        * (m + x_beta + 0.5 * (dm + x_alpha - x_beta))
        * (dn + y_alpha - y_beta)
    )


def magnetic_hamiltonian(
    model: TBModel,
    p: int,
    q: int,
    k1_mag: float,
    k2_mag: float,
    sparse: bool = False,
) -> np.ndarray | sp.csr_matrix:
    """Build ``H_B(k_mag; phi=p/q)`` in lattice-coordinate Landau gauge."""
    if q <= 0:
        raise ValueError("q must be positive")
    phi = p / q
    dim = q * model.n_orb
    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []
    dense = None if sparse else np.zeros((dim, dim), dtype=complex)
    positions = model.positions

    for dm, dn, mat in model.hoppings:
        for mu in range(q):
            mu_raw = mu + dm
            mu_to = mu_raw % q
            cell_crossing = floor(mu_raw / q)
            boundary_phase = np.exp(1.0j * k1_mag * cell_crossing * q)
            k2_phase = np.exp(1.0j * k2_mag * dn)
            for alpha in range(model.n_orb):
                row = mu_to * model.n_orb + alpha
                for beta in range(model.n_orb):
                    amp = mat[alpha, beta]
                    if amp == 0.0:
                        continue
                    phase = np.exp(
                        1.0j
                        * peierls_phase(
                            mu,
                            0,
                            dm,
                            dn,
                            phi,
                            orbital_from=tuple(positions[beta]),
                            orbital_to=tuple(positions[alpha]),
                        )
                    )
                    value = amp * boundary_phase * k2_phase * phase
                    col = mu * model.n_orb + beta
                    if sparse:
                        rows.append(row)
                        cols.append(col)
                        data.append(value)
                    else:
                        dense[row, col] += value

    if sparse:
        out = sp.coo_matrix((data, (rows, cols)), shape=(dim, dim), dtype=complex).tocsr()
        return 0.5 * (out + out.getH())
    assert dense is not None
    return 0.5 * (dense + dense.conj().T)


def magnetic_kmesh(q: int, kmesh: tuple[int, int]) -> np.ndarray:
    """Return a regular magnetic-Brillouin-zone mesh."""
    nk1, nk2 = kmesh
    if nk1 <= 0 or nk2 <= 0:
        raise ValueError("kmesh entries must be positive")
    k1_values = (np.arange(nk1) + 0.5) * (2.0 * np.pi / q) / nk1
    k2_values = (np.arange(nk2) + 0.5) * (2.0 * np.pi) / nk2
    return np.array([(float(k1), float(k2)) for k1 in k1_values for k2 in k2_values], dtype=float)


def spectrum_at_flux(
    model: TBModel,
    p: int,
    q: int,
    kmesh: tuple[int, int] = (1, 12),
    method: str = "dense_full",
    energy_window: tuple[float, float] | None = None,
    n_eigs: int = 64,
    return_eigenvectors: bool = True,
) -> HofstadterSpectrum:
    """Solve magnetic levels at flux ``p/q`` on a magnetic k mesh."""
    k_points = magnetic_kmesh(q, kmesh)
    values: list[np.ndarray] = []
    vectors: list[np.ndarray] = []
    for k1, k2 in k_points:
        if method == "dense_full":
            h = magnetic_hamiltonian(model, p, q, k1, k2, sparse=False)
            if return_eigenvectors:
                evals, evecs = la.eigh(h, check_finite=False)
                vectors.append(evecs)
            else:
                evals = la.eigvalsh(h, check_finite=False)
        elif method == "sparse_window":
            if energy_window is None:
                raise ValueError("energy_window is required for sparse_window")
            center = 0.5 * (energy_window[0] + energy_window[1])
            h_sparse = magnetic_hamiltonian(model, p, q, k1, k2, sparse=True)
            evals, evecs = spla.eigsh(h_sparse, k=n_eigs, sigma=center, which="LM")
            order = np.argsort(evals)
            evals = evals[order]
            evecs = evecs[:, order]
            mask = (evals >= energy_window[0]) & (evals <= energy_window[1])
            evals = evals[mask]
            evecs = evecs[:, mask]
            if return_eigenvectors:
                vectors.append(evecs)
        else:
            raise ValueError("method must be 'dense_full' or 'sparse_window'")
        values.append(np.asarray(evals, dtype=float))

    if method == "dense_full":
        eigenvalues = np.vstack(values)
    else:
        max_len = max(len(v) for v in values)
        eigenvalues = np.full((len(values), max_len), np.nan, dtype=float)
        for i, row in enumerate(values):
            eigenvalues[i, : len(row)] = row
    eigenvectors = np.stack(vectors, axis=0) if return_eigenvectors and method == "dense_full" else None
    return HofstadterSpectrum(
        p=p,
        q=q,
        k_points=k_points,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        state_weight_per_original_cell=1.0 / (q * len(k_points)),
    )


def magnetic_field_tesla_from_q(q: int, flux_quantum_field_tesla: float = 200.0, p: int = 1) -> float:
    """Convert ``p/q`` flux per moire cell to an approximate magnetic field."""
    return flux_quantum_field_tesla * p / q


def rational_flux_grid_for_b_field(
    *,
    b_min: float,
    b_max: float,
    n_points: int,
    q_max: int,
    flux_quantum_field_tesla: float = 200.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Approximate a uniform magnetic-field grid by reduced fluxes ``p/q``.

    The target fields are uniform in ``B``.  Each target is converted to
    ``phi=B/flux_quantum_field_tesla`` and approximated by the nearest reduced
    rational with denominator at most ``q_max``.  Duplicate rationals are
    removed, and the result is sorted by actual magnetic field.
    """
    if b_min <= 0.0 or b_max <= 0.0:
        raise ValueError("b_min and b_max must be positive")
    if b_min >= b_max:
        raise ValueError("b_min must be smaller than b_max")
    if n_points < 2:
        raise ValueError("n_points must be at least 2")
    if q_max < 1:
        raise ValueError("q_max must be positive")
    if flux_quantum_field_tesla <= 0.0:
        raise ValueError("flux_quantum_field_tesla must be positive")

    pairs: dict[tuple[int, int], float] = {}
    for target_b in np.linspace(float(b_min), float(b_max), int(n_points)):
        target_phi = float(target_b) / float(flux_quantum_field_tesla)
        frac = Fraction(target_phi).limit_denominator(int(q_max))
        if frac.numerator <= 0:
            frac = Fraction(1, int(q_max))
        pairs[(int(frac.numerator), int(frac.denominator))] = magnetic_field_tesla_from_q(
            int(frac.denominator),
            flux_quantum_field_tesla=float(flux_quantum_field_tesla),
            p=int(frac.numerator),
        )

    ordered = sorted(pairs.items(), key=lambda item: item[1])
    p_values = np.asarray([key[0] for key, _ in ordered], dtype=int)
    q_values = np.asarray([key[1] for key, _ in ordered], dtype=int)
    b_values = np.asarray([value for _, value in ordered], dtype=float)
    return p_values, q_values, b_values
