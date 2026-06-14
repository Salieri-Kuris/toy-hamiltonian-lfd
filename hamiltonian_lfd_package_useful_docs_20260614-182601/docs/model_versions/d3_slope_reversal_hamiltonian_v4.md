# D3 slope-reversal Hamiltonian v4

Date: 2026-06-12

This version is an aggressive continuation of
`d3_strong_bending_hamiltonian_v3.md`.  It keeps the same D3 two-band parent
structure and still only retunes the scalar `d0(k)` harmonics, but it lowers the
particle-side Lifshitz scale and adds an explicit finite-B beta hybridization
gap in the DOS pipeline.

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
u1 = -0.11824074424096692
u2 = -0.03936296247689948
u3 = -0.023118973634458233
u4 = -0.018372676691527402
u5 = 0.010154034407805107
u6 = 0.014014999412566741
u7 = 0.0012894556500397858
internal_degeneracy = 2
copy_degeneracy = 3
```

## Zero-Field Diagnostics

- Representative red-line Dirac roots on `K'-Gamma-K` (`k2=-k1`):
  `(k1,k2)=(-1.899562633007797, 1.899562633007797)` and
  `(1.899562633007797, -1.899562633007797)`.
- Selected Dirac energy: `E_D=-0.01406456719546756`.
- Selected velocities: `v_parallel=0.5008097114561973`,
  `v_perp=0.5952940449807653`.
- Nearest upper saddle: `k=(-0.04020423873511714, 0.08657713191962246)`.
- Nearest upper saddle energy: `0.008365745279426197`.
- Nearest upper saddle offset from Dirac:
  `Delta_s=0.022430312474893757`.
- Nearest upper saddle Hessian determinant:
  `-3.712084416639034e-08`.
- Nearest lower saddle: `k=(-3.1415564340328412, 4.380271426462912e-06)`.
- Nearest lower saddle energy: `-0.33660940934915917`.
- Nearest lower saddle offset below Dirac:
  `E_D-E_lower=0.3225448421536916`.

## LFD Settings

The default `examples/plot_d3_two_band_lfd.py` run now uses
`--candidate slope_reversal` with

```text
Gamma_alpha = 0.003
Gamma_beta = 0.045
w_s = 0.050
```

Particle-side beta reservoir:

```text
A_beta = 1.10
E_beta = E_s
w_beta = 0.085
B_on = 4.2 T
B_on_width = 0.32 T
```

Finite-B hybridization gap in `D(E,B)`:

```text
gap_depth = 0.65
gap_width = 0.030
gap_center_slope = 0.015 energy/T
gap_B_on = 4.3 T
gap_B_on_width = 0.35 T
```

This final gap is phenomenological.  It represents a beta-sector avoided
crossing / magnetic-breakdown notch in the DOS, not a manually drawn density
line.  The full pipeline is still

```text
D(E,B) -> n(mu,B) -> mu(n,B) -> D(n,B)
```

## Current Status

This version can produce particle-side dark-valley tracks whose fitted high-B
slope is negative in the post-processed DP diagnostic.  The extraction is still
not robust enough to use all overlaid white branches as final physical labels:
the main visual result should be judged from the raw DOS map and from targeted
single-branch diagnostics.
