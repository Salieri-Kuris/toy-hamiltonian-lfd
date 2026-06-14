# D3 strong-bending Hamiltonian v3

Date: 2026-06-12

This version follows `strong_bending_landau_fan_hamiltonian_search_manual.md`.
It keeps the same D3 two-band parent structure and leaves the mass/mixing
sector unchanged, but retunes the scalar `d0(k)` harmonics through `u7` to make
the particle side non-isolated while keeping the hole side isolated.

## Hamiltonian

```text
H(k) = d0(k) sigma_0 + M(k) sigma_z + lambda g(k) sigma_x
```

with

```text
c1 = cos(k1) + cos(k2) + cos(k1 + k2)
c2 = cos(2 k1) + cos(2 k2) + cos(2(k1 + k2))
g  = sin(k1 + k2)
M  = m0 + m1 c1 + m2 c2 + m3 c1^2
d0 = u0 + u1 c1 + u2 c2 + u3 c1^2 + u4 c1 c2
     + u5 c2^2 + u6 c1^3 + u7 c1^2 c2
```

## Parameters

```text
m0 = -0.1
m1 = 0.2
m2 = -0.05
m3 = 0.0
lambda = 0.45
u0 = 0.0
u1 = -0.11823260573457668
u2 = -0.0393635379740415
u3 = -0.023103272118251236
u4 = -0.018348877905879905
u5 = 0.010152092446889431
u6 = 0.014070201232323309
u7 = 0.0013351112620801405
internal_degeneracy = 2
copy_degeneracy = 3
```

## Zero-Field Diagnostics

- Representative red-line Dirac roots on `K'-Gamma-K` (`k2=-k1`):
  `(k1,k2)=(-1.899562633007797, 1.899562633007797)` and
  `(1.899562633007797, -1.899562633007797)`.
- Selected Dirac energy: `E_D=-0.01406583966604184`.
- Selected velocities: `v_parallel=0.5008097114561973`,
  `v_perp=0.5952940449807653`.
- Nearest upper saddle: `k=(0.08095771140963315, 0.08124754961111513)`.
- Nearest upper saddle energy: `0.011394672495507019`.
- Nearest upper saddle offset from Dirac:
  `Delta_s=0.02546051216154886`.
- Nearest upper saddle Hessian determinant:
  `-1.6828090960541768e-06`.
- Nearest lower saddle: `k=(-3.1415564340328412, 4.380271426462912e-06)`.
- Nearest lower saddle energy: `-0.3366106818197334`.
- Nearest lower saddle offset below Dirac:
  `E_D-E_lower=0.32254484215369156`.

## LFD Settings

The default D3 LFD plot now uses this candidate with a stronger, explicitly
finite-field beta-window broadening:

```text
Gamma_alpha = 0.003
Gamma_beta = 0.035
w_s = 0.042
```

The remaining dense-beta reservoir is restricted to the particle side and turns
on around the experimental bending field scale:

```text
A_beta = 0.95
E_beta = E_s
w_beta = 0.095
B_star = 4.5 T
B_on = 4.5 T
B_on_width = 0.40 T
```

This keeps the calculation Hamiltonian-first: the very flat upper saddle lowers
the particle-side Lifshitz scale into the target window, while the finite-field
reservoir mimics magnetic-breakdown or remote-subband hybridization that should
not contaminate the hole side.
