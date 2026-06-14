# Four switchable phenomenology modes: route A+B + visible oscillation

Source: `outputs/reverse_designed_hamiltonian/reverse_candidate_C_pilot_highBkmesh_weak_20260614-214908.npz`

Modes:
- pure: route A+B off, visible oscillation off
- AB_only: route A+B on, visible oscillation off
- visible_only: route A+B off, visible oscillation on
- AB_plus_visible: route A+B on, visible oscillation on

route A+B: `{'Gamma0': 0.0035, 'a_sqrtB': 0.00045, 'lambda_D': 0.0045}`

visible oscillation: `{'smooth_E': 0.018, 'smooth_B': 0.8, 'A_min': 0.18, 'alpha': 1.35}`

Density inversion uses `D_density`; plotted Rxx proxy uses `D_vis(mu(n,B),B)`.

Density DOS integral max relerr by mode: `{'pure': 0.0, 'AB_only': 0.0, 'visible_only': 0.0, 'AB_plus_visible': 0.0}`

Figures:
- `figures/reverse_candidate_C_four_mode_visible_oscillation_EB_20260614-232203.png`
- `figures/reverse_candidate_C_four_mode_visible_oscillation_nB_20260614-232203.png`
- `figures/reverse_candidate_C_four_mode_visible_oscillation_diagnostics_20260614-232203.png`

Data: `outputs/reverse_designed_hamiltonian/reverse_candidate_C_four_mode_visible_oscillation_20260614-232203.npz`
