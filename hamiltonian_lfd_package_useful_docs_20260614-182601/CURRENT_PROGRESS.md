# Current Progress: D3 Toy Hamiltonian Landau Fans

Date: 2026-06-14

## Goal

Construct a D3-symmetric toy Hamiltonian with:

- a Dirac cone on the desired symmetry line,
- particle-hole asymmetric Landau-level spacing,
- a hole-like saddle pocket that can bend particle-side fan lines,
- no dominant dense particle-like pocket that produces the wrong visible fan.

## What Was Tried

1. **Strong-bending v3 model**
   - Produced useful reference plots.
   - Had an unwanted particle-side saddle / particle-like pocket.

2. **Candidate A/B tuning**
   - Adjusted scalar terms to try to turn the problematic pocket into a
     hole-like saddle while preserving broad spacing structure.
   - Candidate B added high-order D3 scalar shape terms.
   - The approach was still too close to the previous model and did not cleanly
     solve the pocket ordering.

3. **Reverse-designed candidate C**
   - Constructed from symmetry and desired spectral requirements rather than by
     simply flipping signs or locally perturbing v3.
   - Included a magnetic-breakdown-inspired replica/splitting postprocess:
     `t_max=0.05`, `onset_field=4.8 T`, `onset_width=0.6 T`,
     `copy_degeneracy=3`.
   - Key source output:
     `outputs/reverse_designed_hamiltonian/reverse_candidate_C_20260614-210842.npz`.

4. **New DOS phenomenology**
   - Implemented route A+B broadened-DOS proxy:
     `Gamma(E,B)^2 = Gamma0^2 + (a_sqrtB sqrt(B))^2
     + (lambda_D D0bar(E))^2`.
   - Reduced route A+B broadening after particle-side fan lines became too
     washed out.
   - Added a separate visible oscillation filter with independent switches:
     `pure`, `AB_only`, `visible_only`, `AB_plus_visible`.
   - Density inversion uses `D_density`; the displayed proxy uses
     `D_vis(mu(n,B),B)`, so visible filtering does not remove states from the
     density count.

## Current Diagnosis

The present issue is model-level, not just phenomenological.

The particle-side upper Dirac branch overlaps in energy with several dense
pockets.  Because scalar route B broadening depends only on a DOS profile
`D0bar(E)`, it cannot distinguish an upper Dirac branch from nearby dense
pocket states at the same energy.  Large particle-side DOS therefore broadens
both the unwanted pocket fan and the desired upper Dirac fan.

The visible oscillation filter helps separate smooth background from oscillatory
visibility, but it also shows that the Hamiltonian itself still needs work:
the desired clear particle-side Dirac fan is not robust once the overlapping
pocket DOS is included.

## Most Relevant Figures

- `figures/reverse_candidate_C_newDOS_weaker_existing_density_20260614-214406.png`
  - Clear weakened route A+B density-axis reference from the full candidate C
    source.
- `figures/reverse_candidate_C_fullsource_four_modes_Dvis_nB_20260614-232834.png`
  - Four switchable visible `D(n,B)` modes using the same full candidate C
    source.
- `figures/reverse_candidate_C_fullsource_four_modes_Ddensity_nB_20260614-232834.png`
  - Four switchable density-control `D(n,B)` modes using the same full source.

## Next Direction

Do not keep tuning route A+B broadening as the main fix.  The next Hamiltonian
iteration should reduce the energy overlap between the upper Dirac branch and
the dense particle-side pockets, while preserving:

- the hole-side Dirac black-line structure,
- particle-hole asymmetric spacing,
- the magnetic-breakdown-inspired particle-side spacing behavior,
- a hole-like saddle pocket capable of contributing smooth background bending.

The visible oscillation filter should remain as a diagnostic layer, with route
A+B and visible filtering independently switchable.

