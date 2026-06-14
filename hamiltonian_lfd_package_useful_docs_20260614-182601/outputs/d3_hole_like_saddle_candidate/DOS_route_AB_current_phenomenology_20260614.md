# Current DOS phenomenology: route A+B

Input data:

```text
source = outputs/d3_hole_like_saddle_candidate/d3_candidate_B_shape_term_fixed_lfd_20260614-204516.npz
model = candidate B shape-term fixed
phenomenology in source = off
```

The manual's route A+B is implemented as a broadened-DOS proxy. It does not add
extra beta-reservoir states, hybridization notches, hand-drawn fan lines, or a
separate visibility mask.

## Route A: global and field-dependent broadening

Starting from the Hofstadter DOS grid `D_Hof(E, B)`, each magnetic-field row is
convolved along the energy axis with a normalized Gaussian kernel:

```text
Gamma_A(B)^2 = Gamma0^2 + GammaB(B)^2
GammaB(B) = a_sqrtB * sqrt(B)
```

The convolution kernel is normalized so the integrated number of states is
preserved row by row.

## Route B: scalar energy-dependent broadening

To mimic stronger quantum-lifetime broadening in high-DOS / heavy-pocket
regions, a scalar DOS profile is extracted from the input data:

```text
D0bar(E) = normalized smooth B-average of D_Hof(E, B)
```

The energy-dependent contribution is

```text
Gamma_D(E) = lambda_D * D0bar(E)
```

Numerically, the DOS is split into energy bins. The DOS weight in each bin is
convolved with a Gaussian whose width is evaluated at that bin center:

```text
Gamma_AB(E_bin, B)^2 = Gamma0^2 + (a_sqrtB * sqrt(B))^2
                     + (lambda_D * D0bar(E_bin))^2
```

This follows the manual's prescription that the width is attached to the level
or energy packet being broadened, not to the output energy coordinate. The
purpose is to keep the state count approximately conserved while washing out
fine oscillations in scalar high-DOS regions.

## Density inversion and proxy

The same broadened total DOS is used for density inversion and for the plotted
proxy:

```text
n(E, B) = integral_{E_D}^{E} D_AB(E', B) dE'
mu(n, B) = inverse of n(E, B)
Rxx_proxy(n, B) proportional to D_AB(mu(n, B), B)
```

Thus the background DOS still contributes to `mu(n, B)` and can bend fan-line
positions, but no extra states are inserted by hand.

## Parameter scan

Three route A+B settings were rendered:

```text
mild:
  Gamma0 = 0.006
  a_sqrtB = 0.0008
  lambda_D = 0.010

medium:
  Gamma0 = 0.008
  a_sqrtB = 0.0012
  lambda_D = 0.018

strong:
  Gamma0 = 0.010
  a_sqrtB = 0.0018
  lambda_D = 0.030
```

The recommended first inspection is `medium`, because it damps the dense
particle-side texture without completely washing out the lower-DOS fan
structure.
