# D3 two-band tight-binding model 与 Peierls-Hofstadter Landau fan 操作手册 v5

## 合并说明

本手册将新的 two-band active Hamiltonian 构造指南与通用
Peierls-Hofstadter magnetic-supercell / DOS / density-normalization 计算流程合并。

合并原则：

1. 用 v4 的 two-band same-active-spectrum 构造替换 v3 中旧的三轨道 Dirac+β Hamiltonian 方案。
2. 不再引入显式第三个低能 β band；β 是 active two-band magnetic spectrum 内的 broad/dense sector。
3. 保留 Peierls-Hofstadter builder、磁通归一化、DOS/density normalization、计算优化和 sanity checks。

---


**Project target.** Construct a minimal tight-binding model for a Landau fan whose main observed features are:

1. hole side: fan-line spacing 6;
2. particle side: fan-line spacing 2;
3. particle side: fan lines bend at higher magnetic field;
4. all relevant lines converge near filling
   $$
   \nu_C=-2,
   $$
   i.e. approximately one hole per moiré unit cell.

The new constraint is important: near this filling we should not introduce an explicit third low-energy band just to supply a \(\beta\) orbit. The minimal active model should contain only **two active bands** forming the Dirac cone. Any \(\beta\)-like background must be generated **inside the same active two-band spectrum**, mainly from the particle-side upper branch through Lifshitz / van-Hove / separatrix / magnetic-breakdown physics.

---

## 0. Revision relative to the previous construction

The previous version implicitly used the logic

$$
\text{two bands form Dirac cone}
+\text{third band supplies }\beta\text{ background}.
$$

This is not consistent with the current filling constraint. The revised logic is

$$
\text{two active bands form Dirac cone}
$$

and

$$
\beta
=\text{broad/dense spectral sector inside the active two-band magnetic spectrum}.
$$

Therefore:

- \(\alpha\) and \(\beta\) are **orbit / magnetic-spectrum labels**, not necessarily distinct zero-field band labels.
- The model should not require a third explicit low-energy band in the active filling window.
- Remote bands may still affect the active two-band model indirectly through high-order harmonics and renormalized parameters, but they should be integrated out rather than kept as a third active band.

---

## 1. Minimal active Hilbert space

### 1.1 Single-valley active block

For one valley block, use a two-orbital Bloch Hamiltonian

$$
H_\tau(\mathbf k)
=
d_0(\mathbf k)\sigma_0
+d_z(\mathbf k)\sigma_z
+d_x(\mathbf k)\sigma_x,
$$

where:

- \(\sigma\) acts in the two active orbital / layer / symmetry-sector degrees of freedom;
- \(\tau=\pm\) labels valleys;
- valley \(U(1)\) means there is no term mixing \(\tau=+\) and \(\tau=-\);
- time reversal exchanges valleys:
  $$
  \mathcal T: H_\tau(\mathbf k)\mapsto H_{-\tau}(-\mathbf k).
  $$

A convenient convention is to construct \(H_+(\mathbf k)\), then impose

$$
H_-(\mathbf k)=H_+^*(-\mathbf k).
$$

For the first implementation, compute one valley block and then add the second valley as a separate block or as an appropriate degeneracy weight, depending on whether finite-\(B\) valley splitting is intentionally included.

### 1.2 Internal twofold degeneracy

The observed spacing 6 and 2 require an unresolved internal twofold factor. In the minimal model this may be represented as a spectator flavor \(s=1,2\), so that

$$
H_{\tau,s}(\mathbf k)=H_\tau(\mathbf k).
$$

If this twofold factor is valley-like, then at finite magnetic field time reversal no longer guarantees exact degeneracy at the same \(B\). One must later check that any valley splitting satisfies

$$
\Delta_{\rm valley}(B)\lesssim \Gamma,
$$

where \(\Gamma\) is the relevant disorder / thermal broadening. If the twofold factor is spin-like, check the projected Zeeman / Ising-SOC splitting in the same way.

---

## 2. Momentum coordinates and D3 harmonics

Use triangular moiré reciprocal-lattice coordinates

$$
\mathbf k=(k_1,k_2),
\qquad
k_i\sim k_i+2\pi.
$$

A useful set of lattice harmonics is

$$
c_1(\mathbf k)
=
\cos k_1+\cos k_2+\cos(k_1+k_2),
$$

$$
c_2(\mathbf k)
=
\cos 2k_1+\cos 2k_2+\cos[2(k_1+k_2)],
$$

and

$$
g(\mathbf k)
=
\sin k_1+\sin k_2-\sin(k_1+k_2).
$$

The function \(g\) vanishes on the three symmetry-related lines

$$
k_2=0,
\qquad
k_1=0,
\qquad
k_1+k_2=0.
$$

These are the natural candidate lines for protected Dirac crossings. Higher harmonics may be added, but they should preserve the desired \(D_3\) symmetry unless an explicit symmetry-breaking perturbation is being studied.

---

## 3. Two-band Hamiltonian ansatz

Use

$$
H(\mathbf k)
=
d_0(\mathbf k)\sigma_0
+M(\mathbf k)\sigma_z
+G(\mathbf k)\sigma_x.
$$

The band energies are

$$
E_\pm(\mathbf k)
=
d_0(\mathbf k)
\pm
\sqrt{M(\mathbf k)^2+G(\mathbf k)^2}.
$$

The Dirac cone condition is

$$
M(\mathbf k_D)=0,
\qquad
G(\mathbf k_D)=0.
$$

A minimal choice is

$$
G(\mathbf k)=\lambda g(\mathbf k),
$$

and

$$
M(\mathbf k)
=
m_0+m_1c_1(\mathbf k)+m_2c_2(\mathbf k)+m_3c_1(\mathbf k)^2+
\cdots.
$$

The scalar part is

$$
d_0(\mathbf k)
=
u_1c_1(\mathbf k)+u_2c_2(\mathbf k)+u_3c_1(\mathbf k)^2+u_4c_1(\mathbf k)c_2(\mathbf k)+\cdots.
$$

Here:

- \(M\) and \(G\) mainly control the Dirac crossing and velocities;
- \(d_0\) mainly controls particle-hole asymmetry, Lifshitz structure, and van-Hove placement;
- higher harmonics mimic remote-band renormalization without adding a third active band.

---

## 4. Dirac cone construction on a C2T-invariant line

Take the representative line

$$
k_2=0.
$$

On this line,

$$
G(k_1,0)=0.
$$

Therefore

$$
H(k_1,0)=d_0(k_1,0)\sigma_0+M(k_1,0)\sigma_z,
$$

and the two sectors cross whenever

$$
M(k_1,0)=0.
$$

Near a root \(k_1=k_D\), write

$$
k_1=k_D+q_\parallel,
\qquad
k_2=q_\perp.
$$

Then

$$
H(\mathbf k_D+\mathbf q)
\simeq
E_D\sigma_0
+v_\parallel q_\parallel\sigma_z
+v_\perp q_\perp\sigma_x
+\text{higher-order scalar terms},
$$

where

$$
E_D=d_0(\mathbf k_D),
\qquad
v_\parallel=\partial_{k_1}M(k_1,0)|_{k_1=k_D},
\qquad
v_\perp=\lambda\partial_{k_2}g(k_1,k_2)|_{\mathbf k_D}.
$$

The leading Dirac cone comes from \(M\) and \(G\). The strong particle-hole asymmetry comes from \(d_0\) and higher-order terms in \(M,G\).

### Required check

The crossing should be linear:

$$
|v_\parallel|>0,
\qquad
|v_\perp|>0.
$$

The Dirac Berry phase around a small loop enclosing the crossing should be approximately

$$
\Phi_B\simeq \pi \quad \text{mod }2\pi.
$$

When plotting bands along high-symmetry paths, do not define band identity by local energy sorting. Use wavefunction-overlap tracking. If two states are degenerate or nearly degenerate, track the degenerate subspace projector rather than individual eigenvectors.

---

## 5. Filling constraint: no explicit third active band

The experimental fan origin is near

$$
\nu_C=-2.
$$

This means the active filling window is close to one hole per moiré unit cell. The construction must therefore satisfy:

1. only the two active bands \(E_\pm\) participate near the fan origin;
2. no third explicit band crosses the same energy window;
3. the zero-field filling associated with the Dirac point is calibrated to \(\nu_C=-2\);
4. any additional density of states used to bend the particle-side fan must come from the same active two-band spectrum.

Operationally, define the zero-field filling relative to a reference energy \(E_{\rm ref}\) as

$$
\nu(E)
=
\nu_{\rm offset}
+
\frac{1}{N_k}
\sum_{n=\pm}\sum_{\mathbf k}
\Theta(E-E_n(\mathbf k))\,w_n,
$$

where \(w_n\) includes spin / valley / spectator degeneracy weights as appropriate. Choose \(\nu_{\rm offset}\) so that

$$
\nu(E_D)=\nu_C=-2.
$$

The purpose is not to derive the absolute filling from this minimal model alone; the purpose is to ensure the active two-band model is used consistently around the experimental origin.

---

## 6. Hole-side target: clean unresolved sixfold Dirac-like sector

The lower branch

$$
E_-(\mathbf k)=d_0(\mathbf k)-\sqrt{M(\mathbf k)^2+G(\mathbf k)^2}
$$

should be relatively clean near the Dirac point. The desired zero-field features are:

1. no nearby van-Hove singularity on the hole side;
2. no strong Lifshitz transition close to \(E_D\) on the hole side;
3. no strong same-band separatrix background in the hole-side energy window;
4. three \(C_3\)-related copies and an unresolved internal twofold factor give effective multiplicity
   $$
   g_h=3\times2=6.
   $$

The expected hole-side sequence is

$$
\nu=-3,-9,-15,\ldots,
$$

which corresponds to a sixfold massless-Dirac-like ladder

$$
\nu_h=-6\left(N+\frac12\right).
$$

### Zero-field diagnostics

Compute:

1. constant-energy contours of \(E_-(\mathbf k)\) for energies below \(E_D\);
2. zero-field DOS of \(E_-\);
3. saddle-point locations by solving
   $$
   \nabla_{\mathbf k}E_-(\mathbf k)=0;
   $$
4. Hessian determinant
   $$
   \det \partial_i\partial_j E_-(\mathbf k).
   $$

The hole side should not contain a nearby saddle with

$$
\det \partial_i\partial_j E_-<0
$$

inside the energy range relevant to the observed clean spacing-6 fan.

---

## 7. Particle-side target I: spacing 2 from resolving threefold pocket degeneracy

The upper branch

$$
E_+(\mathbf k)=d_0(\mathbf k)+\sqrt{M(\mathbf k)^2+G(\mathbf k)^2}
$$

should develop three \(C_3\)-related particle-side pockets. Without magnetic breakdown / orbit-network coupling, these pockets would form an unresolved threefold pocket degeneracy, and with internal twofold degeneracy would give

$$
3\times2=6.
$$

The target is that magnetic field resolves the three pocket sectors:

$$
6\rightarrow 2+2+2.
$$

Then a density scan crosses sublevels separated by

$$
2\frac{eB}{h},
$$

which yields the observed odd sequence

$$
\nu=1,3,5,7,\ldots.
$$



### Required magnetic-spectrum check

For particle-side magnetic subbands, identify whether a cluster follows

$$
6\rightarrow2+2+2.
$$

A practical criterion is

$$
|E_{N,\ell}-E_{N,\ell'}|\gg \Gamma,k_BT
$$

for the three split sectors, while each sector remains internally twofold unresolved.

---

## 8. Particle-side target II: same-band beta background for curvature

Spacing 2 alone does not guarantee bending in a true density-field plane. A clean set of resolved twofold LLs still usually gives linear fan trajectories in \(n\)-\(B\), because the density increment is fixed by degeneracy.

Therefore particle-side bending requires extra states between the sharp twofold \(\alpha\)-like levels:

$$
D_{\rm particle}(E,B)
=
D_\alpha(E,B)+D_\beta(E,B).
$$

In the revised two-band model,

$$
D_\beta
$$

must come from the same active two-band spectrum. Possible sources are:

1. a same-band van-Hove singularity on \(E_+\);
2. a same-band Lifshitz transition / neck formation;
3. a separatrix regime where Landau quantization becomes dense;
4. magnetic breakdown between nearby contour pieces in the repeated-zone representation;
5. dense Hofstadter magnetic subbands generated by the upper branch near a saddle.

The density increment between adjacent visible \(\alpha\)-like levels is then

$$
\Delta n_N(B)
=
2\frac{eB}{h}
+
\int_{E_N^\alpha(B)}^{E_{N+1}^\alpha(B)}
D_\beta(E,B)\,dE.
$$

The first term gives spacing 2. The second term bends the fan line.

### Zero-field condition for same-band beta

The upper branch should contain a saddle / neck not too far above the Dirac point:

$$
\nabla_{\mathbf k}E_+(\mathbf k_{\rm vH})=0,
$$

$$
\det \partial_i\partial_jE_+(\mathbf k_{\rm vH})<0.
$$

Its energy should lie in the particle-side window where curvature is observed:

$$
E_{\rm vH,+}-E_D>0,
$$

but should not be so far away that the relevant Landau orbit never samples it in the experimental field range.

A useful dimensionless tuning target is

$$
E_{\rm vH,+}-E_D
\sim
\text{a few particle-side LL spacings at }B=5\text{--}10\,\mathrm T.
$$

This is not a strict equation; it is a tuning criterion.

---

## 9. Parameter tuning strategy

### 9.1 Dirac position and velocity

Tune \(m_0,m_1,m_2,m_3\) so that \(M(k_1,0)=0\) has a target root \(k_D\) on the representative line.

Tune \(\lambda\) and the slope of \(M\) so that

$$
|v_\parallel|,
\quad
|v_\perp|
$$

are not too small. If velocities are too small, the cone becomes too flat and hole-side LLs may become dense or distorted. If velocities are too large, particle-side curvature may only appear at unrealistically high magnetic field.

### 9.2 Particle-hole asymmetry

Tune scalar harmonics \(u_1,u_2,u_3,u_4,\ldots\) to make

$$
E_+(\mathbf k)
$$

strongly non-Dirac-like while keeping

$$
E_-(\mathbf k)
$$

relatively clean.

The practical goal is:

$$
\text{hole side: clean cone-like constant-energy contours},
$$

$$
\text{particle side: three pockets + nearby saddle/neck}.
$$



---

## 10. Zero-field check list

Before doing Peierls-Hofstadter calculation, run the following diagnostics.

### 10.1 Symmetry checks

Numerically verify the target symmetry relations:

$$
H(C_3\mathbf k)=U_{C_3}H(\mathbf k)U_{C_3}^\dagger,
$$

and, for the representative antiunitary line if using a \(C_2\mathcal T\)-real convention,

$$
H(k_1,0)=H(k_1,0)^*.
$$

For valley blocks:

$$
H_-(\mathbf k)=H_+^*(-\mathbf k).
$$

### 10.2 Dirac checks

Locate all solutions of

$$
M(\mathbf k)=0,
\qquad
G(\mathbf k)=0.
$$

For each solution, compute:

1. energy \(E_D\);
2. velocities \(v_\parallel,v_\perp\);
3. Berry phase around a small loop;
4. symmetry copy relation under \(C_3\) and valley time reversal;
5. filling value \(\nu(E_D)\).

Reject models where the relevant Dirac crossing is not tied to

$$
\nu_C=-2.
$$

### 10.3 Fermi-contour checks

Plot constant-energy contours for \(E_-(\mathbf k)\) and \(E_+(\mathbf k)\):

- hole-side energies below \(E_D\);
- particle-side energies above \(E_D\);
- energies near any particle-side saddle.

Required output:

1. hole side: clean Dirac-like contours;
2. particle side: three \(C_3\)-related pockets;
3. particle side: nearby Lifshitz / neck / saddle region;
4. no explicit third band in the same active energy window.

### 10.4 DOS checks

Compute zero-field DOS:

$$
D_0(E)=\frac{1}{N_k}\sum_{\mathbf k,n=\pm}\delta_\eta(E-E_n(\mathbf k)).
$$

Target:

- low DOS near Dirac point;
- no strong hole-side vH peak too close to \(E_D\);
- a particle-side vH / broad DOS feature at experimentally relevant energy.



## 14. References for mechanism, not as mandatory model ingredients

1. L. Onsager, “Interpretation of the de Haas-van Alphen effect,” *Philosophical Magazine* **43**, 1006–1008 (1952).
2. M. H. Cohen and L. M. Falicov, “Magnetic Breakdown in Crystals,” *Physical Review Letters* **7**, 231–233 (1961).
3. P. Moon, Y. Kim, M. Koshino, T. Taniguchi, K. Watanabe, and J. H. Smet, “Nonlinear Landau Fan Diagram for Graphene Electrons Exposed to a Moiré Potential,” *Nano Letters* **24**, 3339–3346 (2024).

The graphene/hBN paper is useful for the mathematical structure “sharp \(\alpha\) levels plus broad/dense \(\beta\) background,” but the present WSe2-oriented minimal model should not copy its explicit multi-miniband origin. Here \(\beta\) must be generated inside the same active two-band spectrum.

---

# Part II. Peierls-Hofstadter 计算流程

以下部分继承自旧手册的 Peierls substitution / Hofstadter 计算说明。这里的输入不再是旧三轨道模型，而是 Part I 中构造出的 two-band zero-field hopping list。


本部分替换旧版的 Peierls-Hofstadter 求解说明。这里的目标不是使用现成软件包黑箱求谱，而是从我们自己构造的 tight-binding hopping list 出发，自己写一个可控的 magnetic supercell / Hofstadter builder，然后把磁场谱输入到 Landau fan pipeline。

核心路线是：

$$
H_0(\mathbf k)
\longrightarrow
\{h_{\Delta m,\Delta n}\}
\longrightarrow
H_B(k_1^{\rm mag},k_2^{\rm mag};\phi=p/q)
\longrightarrow
D(E,B)
\longrightarrow
n(\mu,B)
\longrightarrow
R_{xx}(n,B).
$$

这里“Peierls substitution”是主计算；普通 Onsager quantization 与 magnetic-breakdown network quantization 只作为低场解释、参数诊断和 sanity check。

---

## 5. 磁通归一化与目标磁场范围

每个 moiré primitive unit cell 的无量纲磁通定义为

$$
\phi
=
\frac{\Phi_{\rm cell}}{\Phi_0}
=
\frac{BA_M}{h/e},
$$

其中 $A_M$ 是 moiré unit-cell area。若 triangular moiré lattice 的周期为 $L_M$，则

$$
A_M=\frac{\sqrt3}{2}L_M^2.
$$

对 twisted WSe$_2$，若取 monolayer lattice constant

$$
a\simeq0.33\,\mathrm{nm},
$$

且 twist angle

$$
\theta\simeq3.89^\circ,
$$

则

$$
L_M\simeq \frac{a}{2\sin(\theta/2)}\simeq4.9\,\mathrm{nm},
$$

$$
A_M\simeq20.7\,\mathrm{nm}^2.
$$

因此一个 flux quantum 对应的磁场尺度为

$$
B_\Phi=\frac{\Phi_0}{A_M}
\simeq
\frac{4.135667696\times10^{-15}\,\mathrm{T\,m^2}}
{20.7\times10^{-18}\,\mathrm{m^2}}
\simeq200\,\mathrm{T}.
$$

所以在本模型中可以近似使用

$$
\phi\simeq \frac{B[\mathrm T]}{200}.
$$

### 5.1 目标磁场 $0.1$--$10$ T 对应的 $\phi$

| $B$ (T) | $\phi\simeq B/200$ | 推荐 rational flux | $q$ if $p=1$ | 计算等级 |
|---:|---:|---:|---:|---|
| 0.1 | 0.0005 | $1/2000$ | 2000 | 极重；不建议全谱 dense diagonalization |
| 0.2 | 0.0010 | $1/1000$ | 1000 | 重；优先 sparse / KPM / 少量能窗 |
| 0.4 | 0.0020 | $1/500$ | 500 | 中重；可做 sparse DOS |
| 0.5 | 0.0025 | $1/400$ | 400 | 中重 |
| 0.8 | 0.0040 | $1/250$ | 250 | 可做机制检查 |
| 1.0 | 0.0050 | $1/200$ | 200 | 可做 sparse / 部分 dense |
| 2.0 | 0.0100 | $1/100$ | 100 | 推荐主扫描点 |
| 3.3 | 0.0167 | $1/60$ | 60 | 推荐主扫描点 |
| 5.0 | 0.0250 | $1/40$ | 40 | 推荐主扫描点 |
| 6.7 | 0.0333 | $1/30$ | 30 | 推荐主扫描点 |
| 8.0 | 0.0400 | $1/25$ | 25 | 推荐主扫描点 |
| 10.0 | 0.0500 | $1/20$ | 20 | 推荐主扫描点 |

第一轮机制验证不需要均匀的 $B$ 网格。建议先取

$$
\phi\in
\left\{
\frac1{200},
\frac1{150},
\frac1{120},
\frac1{100},
\frac1{80},
\frac1{60},
\frac1{50},
\frac1{40},
\frac1{30},
\frac1{25},
\frac1{20}
\right\},
$$

对应大约

$$
B\simeq1\text{--}10\,\mathrm T.
$$

随后再补充

$$
q=250,400,500,1000,2000
$$

检查 $0.1$--$1$ T 的低场趋势。低场 $q$ 很大，因此这里不应一开始就要求完整 dense spectrum；应采用 sparse eigensolver、KPM 或半经典外推。

---

## 6. 为什么不直接在 $0.1$ T 做完整 Hofstadter dense diagonalization

若 zero-field single-valley orbital 数为 $N_{\rm orb}$，取 $q$ 倍 magnetic unit cell 后，单 valley magnetic Bloch Hamiltonian 维度为

$$
N_{\rm mag}=qN_{\rm orb}.
$$

若显式保留 twofold internal flavor，则

$$
N_{\rm mag}=2qN_{\rm orb}.
$$

若再显式保留两个 valleys，则再乘以 $2$。

例如本手册最小 active model 有 Dirac-sector 两个 orbitals 加一个 $\beta$ orbital，因此

$$
N_{\rm orb}=3.
$$

则：

| $q$ | single valley dimension $3q$ | with twofold flavor $6q$ | dense diagonalization |
|---:|---:|---:|---|
| 20 | 60 | 120 | 轻松 |
| 100 | 300 | 600 | 可行 |
| 200 | 600 | 1200 | 可行但多 $k$ 点较慢 |
| 500 | 1500 | 3000 | dense 全谱开始昂贵 |
| 1000 | 3000 | 6000 | 不建议 dense 全谱 |
| 2000 | 6000 | 12000 | 只做 sparse / KPM / 能窗 |

因此，$0.1$ T 对应 $q\sim2000$，直接对每个 $k_{\rm mag}$ 做 dense full diagonalization 并扫很多 $B$ 点通常不划算。更合理的策略是：

1. $1$--$10$ T：用 $q\le200$ 做主机制验证；
2. $0.4$--$1$ T：用 $q=200$--$500$ 做 sparse / KPM DOS；
3. $0.1$--$0.4$ T：用 semiclassical / sparse low-energy windows 做趋势检查，而不是完整 Hofstadter butterfly。

---

## 7. 从 Bloch Hamiltonian 到 hopping list

将 single-valley Bloch Hamiltonian 写成有限 Fourier series：

$$
h(\mathbf k)=\sum_{\Delta m,\Delta n}
h_{\Delta m,\Delta n}
\exp[i(k_1\Delta m+k_2\Delta n)].
$$

这里

$$
\mathbf R=\Delta m\,\mathbf a_1+\Delta n\,\mathbf a_2.
$$

对应 real-space Hamiltonian 是

$$
H_0=
\sum_{m,n}
\sum_{\Delta m,\Delta n}
 c_{m+\Delta m,n+\Delta n}^{\dagger}
 h_{\Delta m,\Delta n}
 c_{m,n}.
$$

代码中应把所有 hopping 存成列表：

```python
hoppings = [
    # (dm, dn, matrix)
    (0, 0, H_onsite),
    (1, 0, T_10),
    (0, 1, T_01),
    (1, 1, T_11),
    # ... include Hermitian conjugates or let builder add them
]
```

注意：

1. 每个 `matrix` 是 $N_{\rm orb}\times N_{\rm orb}$ complex matrix；
2. 若列表只存一半 hoppings，builder 必须自动加入 Hermitian conjugate；
3. 对有限 Fourier harmonics，hopping range 是有限的，Peierls substitution 可直接逐 hopping 加相位；
4. 如果显式保留 valley $U(1)$，可以先只构造 single-valley block，再由 valley relation 生成另一个 block。

---

## 8. Lattice-coordinate Landau gauge 与 Peierls phase

取 unit cell 坐标

$$
\mathbf r=m\mathbf a_1+n\mathbf a_2.
$$

选 lattice-coordinate Landau gauge，使 hopping

$$
(m,n)\rightarrow(m+\Delta m,n+\Delta n)
$$

获得 Peierls phase

$$
\varphi_{m,n}^{\Delta m,\Delta n}
=2\pi\phi
\left(m+\frac{\Delta m}{2}\right)
\Delta n.
$$

于是

$$
h_{\Delta m,\Delta n}
\rightarrow
h_{\Delta m,\Delta n}
\exp\left[i\varphi_{m,n}^{\Delta m,\Delta n}\right].
$$

这个规范的 plaquette lattice curl 是

$$
\sum_{\partial \square}\varphi=2\pi\phi.
$$

若 unit cell 内 orbital 位置不同，orbital coordinate 不能忽略。设 orbital $\alpha$ 在 unit cell 内的位置为

$$
\boldsymbol\rho_\alpha=x_\alpha\mathbf a_1+y_\alpha\mathbf a_2.
$$

则 hopping $\beta\to\alpha$ 的 phase 应使用连续坐标：

$$
\varphi_{m,n,\beta\to\alpha}^{\Delta m,\Delta n}
=2\pi\phi
\left[m+x_\beta+\frac{\Delta m+x_\alpha-x_\beta}{2}\right]
(\Delta n+y_\alpha-y_\beta).
$$

第一版可以先把所有 orbitals 放在 moiré cell center，即

$$
x_\alpha=y_\alpha=0,
$$

但最终必须检查 orbital-position convention 对 spectrum 的影响。

---

## 9. Magnetic Bloch Hamiltonian builder

先采用 $p=1$，即

$$
\phi=\frac1q.
$$

将 magnetic unit cell 沿 $\mathbf a_1$ 方向扩大 $q$ 倍。定义 magnetic subcell index

$$
\mu=0,1,\dots,q-1.
$$

每个 magnetic Bloch momentum 下的 basis 是

$$
|\mu,\alpha\rangle,
\qquad
\mu=0,\dots,q-1,
\quad
\alpha=1,
\dots,N_{\rm orb}.
$$

magnetic Brillouin zone 为

$$
k_1^{\rm mag}\in[0,2\pi/q),
\qquad
k_2^{\rm mag}\in[0,2\pi).
$$

实际 builder 对每个 hopping $(\Delta m,\Delta n,h_{\Delta m,\Delta n})$ 做：

1. 对每个 $\mu$，令

   $$
   \mu_{\rm raw}=\mu+\Delta m.
   $$

2. 计算目标 magnetic subcell

   $$
   \mu' = \mu_{\rm raw}\bmod q.
   $$

3. 计算跨越 magnetic cell 的整数

   $$
   s=\left\lfloor\frac{\mu_{\rm raw}}{q}\right\rfloor.
   $$

4. 乘 magnetic-cell Bloch phase

   $$
   e^{ik_1^{\rm mag}sq}.
   $$

5. 乘 $\mathbf a_2$ 方向 Bloch phase

   $$
   e^{ik_2^{\rm mag}\Delta n}.
   $$

6. 乘 Peierls phase

   $$
   \exp\left[i2\pi\phi\left(\mu+\frac{\Delta m}{2}\right)\Delta n\right].
   $$

7. 把矩阵元加到

   $$
   H_B[(\mu',\alpha),(\mu,\beta)].
   $$

如果使用 general rational flux $\phi=p/q$，只需把上述 phase 中的 $\phi$ 换成 $p/q$。第一轮不建议大量使用 $p>1$，因为 $p>1$ 更容易看到 Hofstadter fractal 细节，而不是清楚的 low-field LL cluster。

---

## 10. 自写代码的结构建议

建议把代码分成四个层级。

### 10.1 Model layer

输入 zero-field hopping list：

```python
@dataclass
class TBModel:
    n_orb: int
    hoppings: list[tuple[int, int, np.ndarray]]
    orbital_pos: np.ndarray | None = None  # shape (n_orb, 2), optional
```

功能：

1. `bloch_hamiltonian(k1, k2)`；
2. `check_hermiticity()`；
3. `check_symmetry_C3_C2T()`；
4. `fermi_contours()`；
5. `extract_hopping_list()` if starting from symbolic Fourier terms.

### 10.2 Magnetic builder layer

核心函数：

```python
def magnetic_hamiltonian(model, p, q, k1_mag, k2_mag, sparse=True):
    ...
```

输出：

$$
H_B(k_1^{\rm mag},k_2^{\rm mag};p/q).
$$

必须支持：

1. dense matrix for small $q$；
2. sparse CSR matrix for large $q$；
3. optional block construction for valley / spin spectator degeneracy；
4. optional orbital projection matrices $P_\alpha,P_\beta$。

### 10.3 Spectrum layer

核心函数：

```python
def spectrum_at_flux(model, p, q, kmesh, method="dense", energy_window=None):
    ...
```

三种求解模式：

1. `dense_full`：小 $q$ 全谱；
2. `sparse_window`：大 $q$ 只算目标能窗附近 eigenvalues；
3. `kpm_dos`：超大 $q$ 直接算 DOS，不显式求所有 eigenvalues。

### 10.4 Fan layer

输入 magnetic spectrum 或 DOS，计算：

$$
D(E,B)\rightarrow n(\mu,B)\rightarrow R_{xx}(n,B).
$$

这里可以复用已有的 `landau_fan_rxx.py` 逻辑，但对于 Hofstadter subbands 更自然的是直接对所有 eigenvalues 做 broadened DOS 和 cumulative density。

---

## 11. 计算优化方案

### 11.1 先做 $1$--$10$ T，不从 $0.1$ T 开始

第一轮主扫描：

$$
q\in\{200,150,120,100,80,60,50,40,30,25,20\}.
$$

这对应

$$
B\simeq1\text{--}10\,\mathrm T.
$$

这一段最相关，因为 particle-side curvature 预期出现在中高场。先确认：

1. hole side 是否保持 spacing $6$；
2. particle side 是否出现 $6\to2+2+2$；
3. 是否存在 $D_\beta(E,B)$ background；
4. bending 是否随打开 $\beta$ coupling 出现。

若这一步失败，不要浪费资源去算 $0.1$ T。

### 11.2 $0.1$--$1$ T 用分层策略

低场对应 $q=200$--$2000$。建议：

- $q=200$--$500$：sparse eigensolver 或 KPM；
- $q=500$--$2000$：优先 KPM DOS 或 semiclassical extrapolation；
- 只在少数代表点检查 full sparse spectrum，不做密集 $B$ 网格。

### 11.3 小 $q$ 用 dense，大 $q$ 用 sparse / KPM

经验 cutoff：

$$
N_{\rm mag}=qN_{\rm orb}N_{\rm flavor}.
$$

- 若 $N_{\rm mag}\lesssim 1500$：可以 dense full diagonalization；
- 若 $1500\lesssim N_{\rm mag}\lesssim 6000$：用 sparse shift-invert 或分能窗 `eigsh`；
- 若 $N_{\rm mag}\gtrsim6000$：不要 full diagonalization，直接用 KPM / Chebyshev DOS。

### 11.4 不要一开始使用大 $k_{\rm mag}$ mesh

对于 low-field LL cluster，magnetic subbands 通常很窄。第一轮可使用

$$
N_{k1}\times N_{k2}=1\times 12,
\quad
2\times 12,
\quad
3\times 18
$$

逐步检查。若 DOS 或 fan-line positions 对 $k$ mesh 不敏感，就无需大 mesh。

对于小 $q$ 或 Hofstadter fractal 明显的 regime，再增大 mesh。

### 11.5 缓存 Peierls phase 与稀疏矩阵结构

对固定 $q$ 和 hopping list，矩阵非零结构不随 $k_1,k_2$ 改变。应缓存：

1. row indices；
2. col indices；
3. hopping orbital matrix entries；
4. Peierls phase factors independent of $k$；
5. boundary crossing integer $s$；
6. $\Delta n$ for $k_2$ phase。

对不同 $k$ 点只更新

$$
e^{ik_1^{\rm mag}sq},
\qquad
 e^{ik_2^{\rm mag}\Delta n}.
$$

这样可避免每个 $k$ 点重新构造完整 hopping loops。

### 11.6 利用 block structure

若 valley $U(1)$ 严格成立，不要在主计算中把两个 valley 放入同一个大矩阵。应分别算：

$$
H_{+,B},
\qquad
H_{-,B},
$$

最后在 DOS 层合并。

若 twofold internal flavor 是完全 spectator degeneracy，也不要显式复制矩阵；直接在 DOS weight 上乘以 $2$。

只有当 Zeeman、orbital magnetic moment 或 internal-flavor splitting 是研究对象时，才显式扩大 Hilbert space。

### 11.7 对 $R_{xx}$ map 使用 streaming accumulation

不要保存所有 $B,k,j$ 的 eigenvectors。若只需要 DOS：

1. 对每个 $B$ 和每个 $k$ 点求 eigenvalues；
2. 立即把 eigenvalues 加入 energy histogram / Gaussian broadened DOS；
3. 释放 eigenvalues；
4. 继续下一个 $k$ 点。

若需要 orbital projection，只保存 projection weights：

$$
w_j^\alpha=\langle u_j|P_\alpha|u_j\rangle,
\qquad
w_j^\beta=\langle u_j|P_\beta|u_j\rangle.
$$

### 11.8 Density normalization

对 flux $p/q$，一个 magnetic unit cell 含 $q$ 个原始 unit cells。若对每个 magnetic $k$ 点有 $qN_{\rm orb}$ 个 bands，则 DOS per original moiré cell 应除以 $q$：

$$
D(E,B)=
\frac{1}{qN_k}
\sum_{\mathbf k_{\rm mag},j}
G_\Gamma(E-E_{j\mathbf k}).
$$

同理，累计 density per original moiré cell 为

$$
n(\mu,B)=
\frac{1}{qN_k}
\sum_{\mathbf k_{\rm mag},j}
F_\Gamma(\mu-E_{j\mathbf k})-n_{\rm ref}(B).
$$

若有 spectator degeneracy $g_{\rm sp}$，则乘上 $g_{\rm sp}$。

这个归一化必须正确，否则 fan spacing 会被错误地放大或缩小。

---

## 12. 从 Peierls spectrum 到目标机制判据

### 12.1 Hole side spacing 6

在 hole side，要求 magnetic spectrum 形成 unresolved sixfold LL clusters。若用 dimensionless filling per moiré cell 表示，则相邻主要 minima 间距应为

$$
\Delta \nu=6.
$$

如果在 Peierls spectrum 中看到 hole side 发生明显

$$
6\to2+2+2,
$$

则说明 hole-side magnetic breakdown 或 symmetry splitting 太强，需要调大 hole-side orbit separation、减弱 hole-side coupling，或增大可见 linewidth。

### 12.2 Particle side spacing 2

在 particle side，目标是原本 sixfold cluster 被解析为三支 twofold subclusters：

$$
6\to2+2+2.
$$

判据：

$$
|E_{N,\ell}-E_{N,\ell'}|
\gtrsim
\Gamma,
\qquad
\ell=0,\pm1.
$$

这会在 density scan 中给出

$$
\nu=1,3,5,7,\dots.
$$

### 12.3 Particle-side curvature

仅有

$$
6\to2+2+2
$$

通常只给 spacing $2$，不保证 $n$-$B$ 中弯曲。弯曲需要额外态数参与 density counting：

$$
\Delta n_N(B)=
2\frac{eB}{h}
+
\int_{E_N}^{E_{N+1}}D_\beta(E,B)\,dE.
$$

所以必须做 control：

1. 关闭 $\beta$ band 或令 $V_{D\beta}=0$；
2. 计算 particle-side fan；
3. 打开 $\beta$ coupling；
4. 比较 bending 是否增强；
5. 输出 $D_\alpha(E,B)$ 与 $D_\beta(E,B)$ 的 orbital-projected DOS。

---

## 13. Semi-classical magnetic breakdown 参数的后验提取

虽然主计算用 Peierls-Hofstadter，但可以从零场 tight-binding band 后验提取 magnetic-breakdown scale。

在 near-touching / avoided crossing 附近，拟合局域两带模型：

$$
H_{\rm loc}(\mathbf q)=
\begin{pmatrix}
\epsilon_1(\mathbf q)&\Delta(\mathbf q)\\
\Delta^*(\mathbf q)&\epsilon_2(\mathbf q)
\end{pmatrix}.
$$

局域 gap 为

$$
2\Delta_{\rm MB}.
$$

两条轨道速度为

$$
\mathbf v_1=\nabla_{\mathbf k}E_1,
\qquad
\mathbf v_2=\nabla_{\mathbf k}E_2.
$$

估计

$$
P_{\rm MB}(E,B)=\exp[-B_0(E)/B],
$$

其中

$$
B_0(E)
\sim
\frac{C\Delta_{\rm MB}^2(E)}
{e\hbar|\mathbf v_1(E)\times\mathbf v_2(E)|}.
$$

$C$ 依赖局域两带 convention，因此第一版只把它用于尺度判断。若 Peierls spectrum 中 splitting onset field 满足

$$
B_{\rm onset}\sim B_0,
$$

则 magnetic-breakdown 解释自洽。若完全不匹配，说明 particle-side spacing $2$ 可能来自 Zeeman / orbital moment / symmetry-breaking splitting，而不是 MB。

---

## 14. 推荐最小计算计划

### Stage 1：zero-field 检查

输出：

1. high-symmetry band structure；
2. $C_2\mathcal T$ line spectrum with symmetry labels；
3. zero-field DOS；
4. hole-side / particle-side Fermi contours；
5. candidate MB junctions and estimated $B_0(E)$。

### Stage 2：主 Peierls-Hofstadter 扫描，$1$--$10$ T

取

$$
q=200,150,120,100,80,60,50,40,30,25,20.
$$

先用 small $k$ mesh 和 dense/sparse 混合方法，输出：

1. $D(E,B)$；
2. $D_\alpha(E,B)$；
3. $D_\beta(E,B)$；
4. $R_{xx}(n,B)$；
5. extracted fan minima。

### Stage 3：低场补点，$0.1$--$1$ T

取

$$
q=250,400,500,1000,2000.
$$

优先使用：

1. KPM / Chebyshev DOS；
2. sparse energy-window eigensolver；
3. semiclassical extrapolation for ordinary LL trend。

低场补点只用于确认 fan lines 是否向 $B\to0$ 正确延伸，不作为第一轮调参主依据。

### Stage 4：control calculations

分别计算：

1. no-$\beta$；
2. no $D$-$\beta$ coupling；
3. full model；
4. increased / decreased $\Gamma_\beta$；
5. valley-resolved DOS；
6. spectator twofold vs explicitly split twofold。

要求：

- hole side 保持 spacing $6$；
- particle side 解析出 spacing $2$；
- particle-side bending 随 $\beta$ sector 和 MB junction 出现；
- 去掉 $\beta$ 或 MB coupling 后 bending 减弱。

---

## 15. 代码级 sanity checks

1. **Hermiticity**：每个 $H_B(k)$ 必须满足

   $$
   \|H_B-H_B^\dagger\|<\epsilon.
   $$

2. **Flux check**：绕一个 primitive plaquette 的 Peierls phase 总和必须为

   $$
   2\pi p/q.
   $$

3. **Gauge check**：更换 Landau gauge 后，sorted spectrum / DOS 不变。

4. **Zero-field envelope**：$q\to\infty$ 时，magnetic subband envelope 回到 zero-field DOS。

5. **Normalization check**：每个 original moiré cell 的总态数应为

   $$
   N_{\rm orb}\times g_{\rm spectator}.
   $$

6. **Degeneracy check**：若某个 spectator twofold 没有显式加入 Hamiltonian，DOS 和 density weight 必须乘 $2$。

7. **Control check**：$w_e=w_o=0$ 时，$D_\beta$ 不应造成 artificial bending。

---

## 16. References

1. L. Onsager, “Interpretation of the de Haas-van Alphen effect,” Philosophical Magazine 43, 1006–1008 (1952).
2. I. M. Lifshitz and A. M. Kosevich, “Theory of magnetic susceptibility in metals at low temperatures,” Soviet Physics JETP 2, 636–645 (1956).
3. M. H. Cohen and L. M. Falicov, “Magnetic Breakdown in Crystals,” Physical Review Letters 7, 231–233 (1961).
4. L. M. Falicov and H. Stachowiak, “Magnetic breakdown in metals,” Physical Review 147, 505–515 (1966).
5. D. R. Hofstadter, “Energy levels and wave functions of Bloch electrons in rational and irrational magnetic fields,” Physical Review B 14, 2239–2249 (1976).
6. P. Moon, Y. Kim, M. Koshino, T. Taniguchi, K. Watanabe, and J. H. Smet, “Nonlinear Landau Fan Diagram for Graphene Electrons Exposed to a Moiré Potential,” Nano Letters 24, 3339–3346 (2024), DOI: 10.1021/acs.nanolett.3c04444.

