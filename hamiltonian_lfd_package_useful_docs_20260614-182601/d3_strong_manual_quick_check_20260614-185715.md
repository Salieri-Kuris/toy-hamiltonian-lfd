# D3 strong manual quick check 20260614-185715

This document records the Hamiltonian and plotting parameters used to reproduce the TB manual quick-check figure:

`figures/tb_hofstadter_lfd_manual_quick_strong_20260614-185715.png`

## Model

Candidate name: `strong`

Code entry point:

- `twse2_continuum_dos/d3_two_band.py`
- Function: `strong_bending_candidate_parameters()`

The Hamiltonian is the D3 two-band model

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

From the reproduction summary:

```text
E_D = -0.014065839666041847
selected Dirac point = (-1.899562633007797, 1.8995626330077968)
v_parallel = 0.3541259430603002
v_perp = 0.6363961030540963
nu_offset = -8.0
selected Dirac filling = -2.0
```

Nearest upper saddle:

```text
k = (1.638793446900877e-05, -1.5996787075867047e-05)
energy = 0.011454488026306848
upper saddle - E_D = 0.025520327692348695
hessian_det = -0.00694925080206236
hessian_trace = 1.1327344457377764
```

Nearest lower saddle:

```text
k = (3.14157819293127, 1.1345340131807546e-05)
energy = -0.33661068200350136
E_D - lower saddle = 0.3225448423374595
hessian_det = -0.03183465114943892
hessian_trace = -0.5454735871324168
```

## Hofstadter / Plotting Parameters

Script:

```text
examples/run_tb_hofstadter_lfd_manual_quick.py
```

Command:

```powershell
C:\ProgramData\anaconda3\python.exe examples\run_tb_hofstadter_lfd_manual_quick.py --candidate strong --q-min 24 --q-max 160 --q-step 4 --kmesh2 7 --eta 0.006 --energy-half-width 0.18 --energy-points 1801 --n-min -0.6 --n-max 1.0 --n-points 1601
```

Numerical settings:

```text
p = 1
q_min = 24
q_max = 160
q_step = 4
q_count = 35
B_range = [1.25, 8.333333333333334] T
kmesh = (1, 7)
eta = 0.006
energy_half_width = 0.18
energy_points = 1801
n_min = -0.6
n_max = 1.0
n_points = 1601
reference_energy = E_D
```

The plotted cyan vertical line in `D(E,B)` is the nearest upper saddle offset:

```text
upper saddle - E_D = 0.025520327692348695
```

## Data Pipeline

This reproduction is a pure Peierls-Hofstadter calculation:

```text
H(k) -> Peierls Hofstadter spectrum -> D(E,B) -> n(mu,B) -> D(n,B)
```

The run explicitly does not add:

```text
beta_reservoir_dos
hand-drawn fan lines
phenomenological magnetic-breakdown replica spectrum
```

## Artifacts

```text
figure:
figures/tb_hofstadter_lfd_manual_quick_strong_20260614-185715.png

data:
outputs/tb_hofstadter_lfd_manual_quick/tb_hofstadter_lfd_manual_quick_strong_20260614-185715.npz

summary:
outputs/tb_hofstadter_lfd_manual_quick/tb_hofstadter_lfd_manual_quick_strong_20260614-185715.json
```

## Notes

This candidate is useful as a reference point because it creates a low particle-side upper saddle at about `E_D + 0.02552` in a symmetry-respecting D3 two-band model. It is not the final desired mechanism: the current broader project goal still requires a Hamiltonian whose hole side is clean and straight while particle-side bending appears only at sufficiently high index.
