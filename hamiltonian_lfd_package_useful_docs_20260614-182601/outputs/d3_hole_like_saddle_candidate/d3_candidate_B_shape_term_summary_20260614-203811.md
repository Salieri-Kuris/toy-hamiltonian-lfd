# Candidate B: high-order D3 d0 shape term

Base: candidate A. Added scalar D3 terms to `d0(k)` with coefficients:

```text
v1 c1 c2^2 + v2 c2^3 + v3 c1^4 + v4 c1^3 c2 + v5 c1^2 c2^2 + v6 c1 c2^3
```

Each basis function is shifted by its value at the Dirac point, so the Dirac energy/root is not directly shifted by the added terms.

```text
v = [-0.00074178902, 0.00041376052, 0.00168148303, -0.00128218385, -0.00107085586, 0.00093484527]
E_D = 0.018906377325
Gamma_upper_minus_ED = 0.040921656652
Gamma_upper_Hessian_eigs = [-0.14155142266658594, 0.11390864293470138]
off_Gamma_particle_min_on_k_minus_k = 0.029992594327 at k=0.412687813022
lower_pi_minus_ED = -0.329963122572
upper_pi_minus_ED = 0.490981031511
TB_hopping_count = 169
TB_hermiticity_error = 1.488e-16
```
