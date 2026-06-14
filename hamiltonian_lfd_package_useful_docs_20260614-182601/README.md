# Toy Hamiltonian LFD Workspace

This snapshot contains the D3 toy-Hamiltonian and Peierls-Hofstadter tools used
to study zero-field pockets, magnetic DOS, and density-axis Landau fan diagrams.

The current conclusion is negative but useful: the latest candidate still has a
model-level problem on the particle side.  The upper Dirac branch overlaps in
energy with several dense pockets, so scalar DOS broadening tends to wash out
the visible particle-side fan lines.  The phenomenological visibility layer can
diagnose this, but it does not by itself fix the Hamiltonian.

## Main Code

- `twse2_continuum_dos/d3_two_band.py`: D3 two-band Hamiltonian helpers and
  finite-range tight-binding conversion.
- `twse2_continuum_dos/d3_physical_reservoir.py`: exploratory physical
  reservoir / pocket diagnostics.
- `twse2_continuum_dos/peierls_hofstadter.py`: generic Peierls-substitution
  Hofstadter solver for finite-range tight-binding models.
- `twse2_continuum_dos/dos.py`: DOS broadening, density regridding, magnetic
  breakdown replicas, low-DOS branch tracking, and phenomenology utilities.
- `twse2_continuum_dos/tb_hofstadter_lfd_workflow.py`: small workflow helper
  for pure Hofstadter DOS/LFD runs.
- `examples/`: plotting and reproduction entry points.
- `tests/`: regression tests for DOS and Hofstadter utilities.

## Current Working Notes

- `CURRENT_PROGRESS.md`: current stage summary and next recommended direction.
- `DOS_as_Rxx_proxy_manual.md`: manual for treating DOS as an `Rxx` proxy,
  including route A+B broadening and visible oscillation filtering.
- `TB_Hofstadter_LFD_mask_manual.md`: Hofstadter/LFD workflow notes.
- `d3_two_band_peierls_landau_fan_manual_v5.md`: two-band D3 Peierls/LFD notes.
- `d3_strong_bending_hamiltonian_v3.md`: older strong-bending v3 model note.

## Important Generated Outputs

The most relevant outputs for the latest discussion are under
`outputs/reverse_designed_hamiltonian/` and `figures/`:

- `reverse_candidate_C_20260614-210842.npz`: full candidate C source data.
- `reverse_candidate_C_newDOS_weaker_from_existing_20260614-214406.npz`:
  weakened route A+B broadening on the full candidate C source.
- `reverse_candidate_C_fullsource_four_mode_visible_oscillation_20260614-232834.npz`:
  four switchable modes on the same full candidate C source:
  `pure`, `AB_only`, `visible_only`, and `AB_plus_visible`.
- `reverse_candidate_C_fullsource_four_modes_Dvis_nB_20260614-232834.png`:
  four-mode visible `D(n,B)` comparison using the full candidate C source.
- `reverse_candidate_C_newDOS_weaker_existing_density_20260614-214406.png`:
  clear reference `D(n,B)` figure from weakened route A+B on the full source.

## Environment

On this computer, use the local Conda base Python:

```powershell
C:\ProgramData\anaconda3\python.exe -m pytest
```

For a generic Python environment:

```bash
python -m pip install -e .
python -m pytest
```

