# TB--Peierls--Hofstadter--LFD 手册：构造零场哈密顿量、检测 Hofstadter spectrum、定义 $\alpha/\beta$ mask

## 0. 目的与边界

本手册的目标是给出一个一致的计算流程：

$$
H_0(\mathbf k)
\longrightarrow
H_B
\longrightarrow
D(E,B)
\longrightarrow
n(\mu,B)
\longrightarrow
\mathrm{LFD}.
$$

这里 $H_0(\mathbf k)$ 是零场平移不变 tight-binding / Bloch Hamiltonian；$H_B$ 是通过 Peierls substitution 得到的 Hofstadter Hamiltonian；$D(E,B)$ 是从同一个 $H_B$ 计算得到的总 DOS；Landau fan diagram 应从总谱、总 DOS、累计 density 推出，而不是预先手动写 fan line。

本手册重点处理三个问题：

1. 如何构造具有 $D_3+\mathcal T+U(1)_{\rm valley}$ 对称性、且 Dirac cone 位于 $C_2$-invariant line 上的零场 TB Hamiltonian；
2. 在绘制 LFD 前，应如何检测零场 Fermi contour、zero-field DOS、Hofstadter spectrum 是否具备产生 particle-side bending 的必要结构；
3. 在完整 Peierls--Hofstadter 框架下，如何定义用于诊断的 $\alpha/\beta$ phase-space mask，以及这种 mask 在什么条件下是 well-defined。

需要明确：本手册不再使用“手动 inter-patch coupling”作为零场 Hamiltonian 的一部分。若某个 effective inter-patch coupling 出现，它只能是从零场 band geometry 与有限磁场 magnetic breakdown 推导出的 semiclassical effective object，而不是 microscopic zero-field hopping term。零场 microscopic Hamiltonian 必须保持原始平移不变性。

---

## 1. 当前物理目标

实验/现象目标可以抽象为：

1. hole side 有较大的 unresolved degeneracy，给近似线性 fan；
2. particle side 的可见 fan spacing 比 hole side 小；
3. particle side fan line 在高场出现强烈弯曲，甚至局部斜率接近变号；
4. 这种强弯曲不应通过手动添加 background DOS 得到，而必须来自同一个 Peierls--Hofstadter spectrum。

在当前模型中，由于 fermion doubling 和 valley $U(1)$，总 Dirac cone 数目会比之前的简化 counting 翻倍。例如，一个 valley 中若有两组 $C_3$-related Dirac roots，则

$$
N_{\rm cone}^{\rm per\ valley}=2\times3=6,
$$

两个 valley 且 valley $U(1)$ 保持时，总共

$$
N_{\rm cone}^{\rm total}=12.
$$

这个绝对数目会影响最终 filling label 的整体归一化。但是本手册的核心不是绝对简并数，而是以下结构：

$$
\text{hole side: large unresolved manifold},
$$

$$
\text{particle side: clean }\alpha\text{ ladder}+\text{dense }\beta\text{ reservoir}.
$$

也就是说，目标是让 particle side 的可见 LL 来自一个较小 degeneracy 的 clean closed orbit，而其他同能区态进入 dense / open / Hofstadter-background-like sector，从而改变从 energy 到 density 的映射。

---

## 2. 对称性与单 valley Bloch Hamiltonian

### 2.1 Valley $U(1)$ 与 time reversal

取 valley label

$$
\tau=\pm.
$$

valley $U(1)$ 意味着零场 Hamiltonian block diagonal：

$$
H_0(\mathbf k)=
\begin{pmatrix}
H_+(\mathbf k)&0\\
0&H_-(\mathbf k)
\end{pmatrix}.
$$

Time reversal 交换两个 valley：

$$
H_-(\mathbf k)=H_+^*(-\mathbf k).
$$

因此可以先构造一个 valley 的 $H_+(\mathbf k)$，另一个 valley 由 $\mathcal T$ 得到。

### 2.2 Representative two-band form

一个最小 two-band valley block 可写为

$$
H_+(\mathbf k)
=
d_0(\mathbf k)\sigma_0
+d_z(\mathbf k)\sigma_z
+d_x(\mathbf k)\sigma_x.
$$

这里 $\sigma_i$ 是 two-band / orbital / $C_2$-sector pseudospin。若希望用 $\sigma_y$ 写 off-diagonal 项，也可以作 unitary rotation；本手册采用 real form，便于表示 $C_2\mathcal T$ 结构。

### 2.3 $C_3$ symmetry

在 triangular moiré reciprocal coordinates $(k_1,k_2)$ 中，可以取

$$
C_3:(k_1,k_2)\mapsto(k_2,-k_1-k_2).
$$

若内部表示先取最小形式

$$
U_{C_3}=\sigma_0,
$$

则要求

$$
d_a(C_3\mathbf k)=d_a(\mathbf k),
\qquad a=0,x,z.
$$

更一般地，若 $U_{C_3}\ne\sigma_0$，则要求

$$
U_{C_3}H_+(\mathbf k)U_{C_3}^\dagger
=H_+(C_3\mathbf k).
$$

### 2.4 $C_2$-invariant line 与 protected crossing

取 representative $C_2$ / mirror-like operation：

$$
C_2:(k_1,k_2)\mapsto(k_1+k_2,-k_2).
$$

其 fixed line 是

$$
k_2=0.
$$

取 internal representation

$$
U_{C_2}=\sigma_z.
$$

要求

$$
d_0(C_2\mathbf k)=d_0(\mathbf k),
$$

$$
d_z(C_2\mathbf k)=d_z(\mathbf k),
$$

$$
d_x(C_2\mathbf k)=-d_x(\mathbf k).
$$

于是沿 fixed line $k_2=0$，有

$$
d_x(k_1,0)=0.
$$

Hamiltonian 在该线上为

$$
H_+(k_1,0)=d_0(k_1,0)\sigma_0+d_z(k_1,0)\sigma_z,
$$

不同 $C_2$ eigenvalue sector 不 hybridize。当

$$
d_z(k_D,0)=0
$$

时产生 symmetry-protected crossing。其 $C_3$ copies 位于另外两条 $C_2$-related lines 上。

---

## 3. 一个可用的 triangular-lattice harmonic 构造

定义

$$
c_1(\mathbf k)=\cos k_1+\cos k_2+\cos(k_1+k_2),
$$

$$
g(\mathbf k)=\sin k_1+\sin k_2-\sin(k_1+k_2).
$$

它们满足：

$$
c_1(C_3\mathbf k)=c_1(\mathbf k),
\qquad
c_1(C_2\mathbf k)=c_1(\mathbf k),
$$

$$
g(C_3\mathbf k)=g(\mathbf k),
\qquad
g(C_2\mathbf k)=-g(\mathbf k).
$$

因此可以取

$$
d_x(\mathbf k)=\lambda g(\mathbf k).
$$

沿 $k_2=0$，

$$
g(k_1,0)=0.
$$

### 3.1 Mass term and fermion doubling

取

$$
d_z(\mathbf k)
=\Delta[c_1(\mathbf k)-c_A][c_1(\mathbf k)-c_B].
$$

在 $k_2=0$ 上，

$$
c_1(k_1,0)=1+2\cos k_1.
$$

如果

$$
-1<c_A,c_B<3,
$$

则在一个 representative half-line 上有两个 Dirac roots：

$$
1+2\cos k_{D,A}=c_A,
$$

$$
1+2\cos k_{D,B}=c_B.
$$

经 $C_3$ 复制后，一个 valley 中有

$$
2\times3=6
$$

个 Dirac cones；两个 valleys 且 valley $U(1)$ 保持时，总共

$$
12
$$

个 Dirac cones。这个 doubling 是模型结构的一部分，但后续重点是 particle/hole side 的可见 spacing 与 bending，而不是绝对 cone number。

### 3.2 Scalar dispersion and particle-hole asymmetry

最小 scalar term 可取

$$
d_0^{(1)}(\mathbf k)=\eta(c_m-c_1(\mathbf k)),
\qquad
c_m=\frac{c_A+c_B}{2}.
$$

要产生强 particle-hole asymmetry，需要加入更高阶 $C_3$-invariant harmonics，例如

$$
c_2(\mathbf k)=\cos(2k_1)+\cos(2k_2)+\cos[2(k_1+k_2)],
$$

$$
c_3(\mathbf k)=\cos(3k_1)+\cos(3k_2)+\cos[3(k_1+k_2)].
$$

一般写成

$$
d_0(\mathbf k)
=d_0^{(1)}(\mathbf k)
+a_2c_2(\mathbf k)
+a_3c_3(\mathbf k)
+a_{11}c_1(\mathbf k)^2+
\cdots.
$$

但如果希望所有 Dirac roots 仍然共享同一 reference energy，应对 added harmonics 加 constraints：

$$
\delta d_0(\mathbf Q_a)=0,
$$

必要时还要求

$$
\nabla\delta d_0(\mathbf Q_a)=0.
$$

这样高阶项主要改变 finite-energy particle/hole contours，而不移动 Dirac crossing energy 或过度改变 Dirac cone 的线性结构。

---

## 4. 产生不同 particle/hole spacing 的正确零场目标

在不手动加入 inter-patch coupling 的 microscopic TB 中，particle/hole spacing difference 必须来自零场 Fermi-surface topology 和有限磁场谱的真实结构。

### 4.1 Hole side target

对 $E<E_D$，目标是每个 valley 中有三组 $C_3$-related disconnected pockets：

$$
\mathcal C_h(E)=
\mathcal C_{h,0}(E)
\cup
\mathcal C_{h,1}(E)
\cup
\mathcal C_{h,2}(E),
$$

with

$$
\mathcal C_{h,a+1}=C_3\mathcal C_{h,a}.
$$

因此 hole side 的 visible degeneracy 是

$$
g_h=3\times g_{\rm valley/internal}.
$$

若 valley 或 fermion doubling 使绝对简并翻倍，只需整体乘上对应 factor；检测重点是 hole side 是否保持 full unresolved manifold。

### 4.2 Particle side target

对 $E>E_D$，如果仍然是每个 valley 三组 disconnected pockets，则 noninteracting $C_3$-preserving 模型会继续给三重 star degeneracy；particle side 不会自然变成 smaller spacing。

因此 particle side 的零场目标应是：

$$
\mathcal C_e(E)=\mathcal C_\alpha(E)\cup\mathcal C_\beta(E),
$$

其中

$$
\mathcal C_\alpha(E)=\text{one }C_3\text{-invariant clean closed orbit per valley},
$$

而

$$
\mathcal C_\beta(E)=\text{large / saddle-connected / near-open / magnetic-breakdown network}.
$$

于是可见 clean LL degeneracy 变为

$$
g_\alpha=1\times g_{\rm valley/internal},
$$

而 $\beta$ sector 不是另一套 clean fan，而是 dense background / reservoir。

这是一条不破坏 $C_3$ 的路线：不是选择 $M_1$ 做 $\alpha$，而是让 particle side 的 closed orbit 本身成为整体 $C_3$-invariant orbit。

---

## 5. 产生 bending 的零场 band-geometry 条件

普通 isolated LL ladder 即使 $E_N(B)$ 非线性，在 $n$-$B$ 平面仍倾向于给 linear fan，因为每个 LL 的态数由 flux degeneracy 决定。因此 bending 需要额外累计态数：

$$
\nu_i(B)=\nu_C+g_\alpha i\phi(B)+\delta\nu_{\beta,i}(B),
$$

where

$$
\delta\nu_{\beta,i}(B)
=\int_{E_{\rm ref}}^{E_i^\alpha(B)}D_\beta(E,B)dE.
$$

要强 bending，需要

$$
\left|\frac{d\delta\nu_{\beta,i}}{dB}\right|
\sim
g_\alpha i\frac{d\phi}{dB}.
$$

因此零场 Hamiltonian 应产生：

1. clear $\alpha$ closed orbit；
2. large phase-space-volume $\beta$ network；
3. particle-side vHS / Lifshitz neck / broad saddle near the $\alpha$ LL energy window；
4. hole side 没有对应强 $\beta$ reservoir。

---

## 6. 绘制 Hofstadter spectrum / LFD 前的零场检测

### 6.1 Symmetry check

检查数值 Hamiltonian 是否满足：

$$
\|H_+(C_3\mathbf k)-U_{C_3}H_+(\mathbf k)U_{C_3}^\dagger\|<\epsilon,
$$

$$
\|H_+(C_2\mathbf k)-U_{C_2}H_+(\mathbf k)U_{C_2}^\dagger\|<\epsilon,
$$

以及

$$
H_-(\mathbf k)=H_+^*(-\mathbf k).
$$

Valley $U(1)$ 要求

$$
H_{+-}=H_{-+}=0.
$$

### 6.2 Dirac cone and $C_2$ label check

沿 representative $C_2$ line，例如 $k_2=0$，检查：

$$
[H_+(k_1,0),U_{C_2}]=0.
$$

对本征态计算

$$
\xi_n(k_1)=\langle u_n(k_1,0)|U_{C_2}|u_n(k_1,0)\rangle.
$$

一个 protected crossing 应表现为两条具有不同 $C_2$ label 的 bands crossing。注意：不能仅使用 energy sorting 判断 band identity。沿路径需要 continuity-based band tracking：

$$
O_{mn}(k_i,k_{i+1})=|\langle u_m(k_i)|u_n(k_{i+1})\rangle|^2,
$$

并用最大 overlap 追踪 band connectivity。在 near-degeneracy 处，应追踪 projector/subspace，而不是单一本征矢。

### 6.3 Hole-side Fermi contour check

对一系列

$$
E=E_D-\delta E
$$

画出

$$
E_n(\mathbf k)=E.
$$

要求：每个 valley 中为三组 disconnected $C_3$-related pockets。

记录：

$$
N_h^{\rm comp}(E)=3
$$

在目标 hole energy window 内稳定成立。

### 6.4 Particle-side Fermi contour check

对一系列

$$
E=E_D+\delta E
$$

画出

$$
E_n(\mathbf k)=E.
$$

目标是：

1. 存在一个 single $C_3$-invariant closed $\alpha$ orbit；
2. 同一能区存在 $\beta$ saddle / large contour / connected network；
3. $\beta$ 不是远离 $\alpha$ LL 能区的 remote feature。

如果 particle side 仍然只是三组 disconnected pockets，则该 noninteracting $C_3$-preserving model 仍会给 star degeneracy，不会自然产生 smaller visible spacing。

### 6.5 Zero-field DOS and reservoir capacity

计算

$$
D_0(E)=\frac{1}{N_k}\sum_{n,\mathbf k}G_\eta(E-E_n(\mathbf k)).
$$

并用后文定义的 masks 估计

$$
D_{0,\alpha}(E),\qquad D_{0,\beta}(E).
$$

要求 particle side target window 内：

$$
D_{0,\beta}(E)\text{ large or peaked},
$$

并且累计容量

$$
C_\beta(E_1,E_2)=\int_{E_1}^{E_2}D_{0,\beta}(E)dE
$$

达到实验所需 bending 的 filling 量级。若目标弯曲相对 linear extrapolation 偏移为 $O(0.2-1)$ moiré filling，则 $C_\beta$ 也必须为同量级。

### 6.6 Magnetic-breakdown scale estimate

在 $\alpha$ orbit 与 $\beta$ network 的 near-touching / avoided-crossing region，估算局部 gap

$$
\Delta_{\rm MB}
$$

与 velocities

$$
\mathbf v_1=\nabla_\mathbf k E_1,
\qquad
\mathbf v_2=\nabla_\mathbf k E_2.
$$

粗略 magnetic-breakdown scale 满足

$$
B_0\propto \frac{\Delta_{\rm MB}^2}{|\mathbf v_1\times\mathbf v_2|}.
$$

目标是

$$
B_0\sim 5\text{--}10\ {\rm T}
$$

或对应实验强 bending 出现的场区。若 $B_0$ 太大，$\beta$ reservoir 在实验场区打不开；若太小，低场 spacing 可能已经被破坏。

---

## 7. Peierls--Hofstadter spectrum 的检测

你已有 Hofstadter / Peierls 代码时，建议先不要直接进入 $D(n,B)$ 或 $R_{xx}$。先看 energy-field DOS：

$$
D(E,B)=\frac{1}{N_k}\sum_{\lambda,\mathbf k_{\rm MBZ}}
G_\Gamma(E-E_{\lambda\mathbf k}(B)).
$$

目标图像是：

$$
D(E,B)=\text{clear }\alpha\text{-LL ridges}+	ext{dense }\beta\text{-background}.
$$

若图像只有稀疏清晰 LL ridges 和黑背景，那么模型仍处在 ordinary LL / Hofstadter regime，通常不会给强 bending。

若 $\beta$ sector 也表现为清晰、分离的 LL ridges，那么它会产生额外 linear fan / gap labels，而不是充当 reservoir。

---

## 8. $\alpha/\beta$ mask 的定义

### 8.1 为什么需要 mask

在 full two-band $H_0(\mathbf k)$ 中，$\alpha$ 与 $\beta$ 通常不是不同 microscopic orbital，而是不同 zero-field phase-space regions / semiclassical orbit families。因此 $D_\alpha,D_\beta$ 一般不是严格的 internal-projected DOS。

正确态度是：

1. 最终 LFD 必须使用 total DOS $D(E,B)$；
2. $D_\alpha,D_\beta$ 只作为诊断工具，用来判断 bending 是否来自 clean $\alpha$ ridge 和 dense $\beta$ reservoir；
3. 如果没有清楚的 zero-field phase-space separation，则 $\alpha/\beta$ mask 不 well-defined。

### 8.2 Zero-field Bloch projector mask

零场对角化：

$$
H_0(\mathbf k)|u_{n\mathbf k}\rangle
=E_n(\mathbf k)|u_{n\mathbf k}\rangle.
$$

在 active energy window 内定义 smooth masks：

$$
\chi_\alpha(n,\mathbf k),\qquad \chi_\beta(n,\mathbf k),
$$

其中

$$
0\le \chi_\eta\le1,
\qquad \eta=\alpha,\beta.
$$

理想情况下

$$
\chi_\alpha+\chi_\beta\simeq1
$$

在 active particle-side window 内；若存在其他无关态，则允许

$$
\chi_\alpha+\chi_\beta<1.
$$

对应的零场 phase-space projectors 是

$$
P_\eta^{(0)}
=
\sum_n\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}
\chi_\eta(n,\mathbf k)
|u_{n\mathbf k}\rangle\langle u_{n\mathbf k}|,
\qquad
\eta=\alpha,\beta.
$$

若存在 near-degenerate band group，必须用 subspace projector：

$$
P_{\mathcal A}(\mathbf k)=\sum_{n\in\mathcal A}|u_{n\mathbf k}\rangle\langle u_{n\mathbf k}|,
$$

而不是单条能带的本征矢。

### 8.3 Mask by contour components

推荐流程：

1. 选定 particle-side active energy slices：
   $$
   E_j\in [E_{p,\min},E_{p,\max}].
   $$
2. 对每个 $E_j$，找等能线 components：
   $$
   E_n(\mathbf k)=E_j.
   $$
3. 标记 clean single $C_3$-invariant closed component 为 $\alpha$。
4. 标记 saddle-connected / large / near-open / high-DOS network components 为 $\beta$。
5. 用 connected-component overlap 在相邻 $E_j$ 之间追踪 component identity。
6. 将 contour labels 延拓到附近 phase-space region，得到 smooth $\chi_\alpha,\chi_\beta$。

一个实用的 smooth mask 形式是：

$$
\tilde\chi_\eta(n,\mathbf k)
=
W_E(E_n(\mathbf k))
\exp\left[-\frac{d_\eta(n,\mathbf k)^2}{2\sigma_k^2}\right],
$$

其中 $d_\eta$ 是 $(n,\mathbf k)$ 到对应 contour family 的最小距离，$W_E$ 是 active energy window function。

归一化为

$$
\chi_\eta=
\frac{\tilde\chi_\eta}{\tilde\chi_\alpha+\tilde\chi_\beta+\epsilon},
\qquad \eta=\alpha,\beta,
$$

或在只关心 $\alpha$ 的情况下取

$$
\chi_\beta=W_E-\chi_\alpha.
$$

---

## 9. 在 Hofstadter eigenstates 上计算 projected DOS

### 9.1 Momentum-space Hofstadter basis

如果你的 Hofstadter Hamiltonian 使用 momentum-space Harper basis，通常 basis component 可标记为

$$
|\kappa_m,\mu\rangle,
$$

其中 $\mu$ 是 orbital / two-band index，$\kappa_m$ 是由 magnetic coupling 连接的一组 crystal momenta。

零场 mask projector 在这个 basis 中近似为 block diagonal：

$$
[P_\eta^{(0)}]_{m\mu,n\nu}
=\delta_{mn}
[P_\eta(\kappa_m)]_{\mu\nu},
$$

其中

$$
P_\eta(\mathbf k)
=
\sum_n\chi_\eta(n,\mathbf k)
|u_{n\mathbf k}\rangle\langle u_{n\mathbf k}|.
$$

对 Hofstadter eigenstate

$$
|\psi_\lambda(B)\rangle
=\sum_{m,\mu}c_{m\mu}^{\lambda}|\kappa_m,\mu\rangle,
$$

计算

$$
w_\eta^\lambda(B)
=\langle\psi_\lambda|P_\eta^{(0)}|\psi_\lambda\rangle.
$$

然后

$$
D_\eta(E,B)
=\sum_\lambda w_\eta^\lambda(B)
G_\Gamma(E-E_\lambda(B)).
$$

### 9.2 Real-space finite torus basis

如果你的 Hofstadter Hamiltonian 使用 real-space basis，可先构造 real-space kernel：

$$
P_\eta^{(0)}(\mathbf R\mu,\mathbf R'\nu)
=
\frac{1}{N_k}
\sum_{\mathbf k,n}
\chi_\eta(n,\mathbf k)
 e^{i\mathbf k\cdot(\mathbf R-\mathbf R')}
 u_{n\mu}(\mathbf k)u^*_{n\nu}(\mathbf k).
$$

然后对 finite-$B$ eigenstate $\psi_\lambda$ 计算

$$
w_\eta^\lambda
=\psi_\lambda^\dagger P_\eta^{(0)}\psi_\lambda.
$$

### 9.3 Projected density

同时定义 projected cumulative density：

$$
N_\eta(\mu,B)
=\sum_\lambda w_\eta^\lambda(B)
F_\Gamma(\mu-E_\lambda(B)),
$$

其中 $F_\Gamma$ 是 Gaussian CDF。该量用于检测 $\beta$ reservoir 是否足够改变 fan-line filling。

---

## 10. Mask 什么时候 well-defined

$\alpha/\beta$ mask 是诊断工具。它 well-defined 需要以下条件。

### 10.1 Zero-field phase-space separation

在目标 energy window 内，$\alpha$ 与 $\beta$ 对应不同 contour families：

$$
\mathcal C_\alpha(E)\cap\mathcal C_\beta(E)
$$

只在少数 near-breakdown / saddle regions 接近，而不是大面积混合。

### 10.2 Topology robust under energy variation

在 $[E_{p,\min},E_{p,\max}]$ 内，$\alpha$ contour identity 可以通过 connected-component tracking 连续追踪。若 Lifshitz transition 太密集，component identity 频繁变化，则 mask 不稳定。

### 10.3 Robust under mask smoothing

改变 $\sigma_k,\sigma_E$ 后，主要诊断量应稳定，例如：

$$
D_\alpha\text{ 是否有 clear ridges},
$$

$$
D_\beta\text{ 是否 dense},
$$

$$
\Delta\nu_{\beta,i}(B)
$$

的量级和趋势不应发生定性改变。

### 10.4 Compatible with magnetic uncertainty

有限磁场下 wavepacket 的 momentum uncertainty 约为

$$
\delta k\sim \ell_B^{-1}.
$$

mask 边界不能比这个尺度更尖锐地决定物理结论。若 $\alpha$ 与 $\beta$ contour 的分离小于 $O(\ell_B^{-1})$，单个 eigenstate 的 $\alpha/\beta$ 归属会变模糊；此时 projected DOS 仍可作为 smooth diagnostic，但不能把单条 Hofstadter level 严格归属为 $\alpha$ 或 $\beta$。

### 10.5 Projected weights should polarize the expected features

成功情形应看到：

$$
w_\alpha^\lambda\approx1
$$

on clean $\alpha$ LL ridges, while dense background has large

$$
w_\beta^\lambda.
$$

若所有 relevant states 都有

$$
w_\alpha^\lambda\sim w_\beta^\lambda\sim\frac12,
$$

则 $\alpha/\beta$ distinction 已不再清楚。

### 10.6 Not an exact spectral decomposition unless internal projectors exist

若 $\alpha,\beta$ 是真实 Hilbert-space subspaces，则 $P_\alpha,P_\beta$ 给严格 projected DOS。若它们只是 zero-field phase-space masks，则分解依赖 mask；最终物理 LFD 必须使用 total DOS。

---

## 11. Hofstadter spectrum 的 bending 检测量

### 11.1 Track $\alpha$ ridges

从 $D_\alpha(E,B)$ 中追踪 clear $\alpha$ LL ridge energies：

$$
E_i^\alpha(B).
$$

不要只按能量排序；应按 ridge continuity / spectral weight continuity 追踪。

### 11.2 Compute beta capacity between alpha ridges

定义

$$
\Delta\nu_{\beta,i}(B)
=
N_\beta(E_{i+1}^\alpha(B),B)-N_\beta(E_i^\alpha(B),B).
$$

若

$$
\Delta\nu_{\beta,i}(B)\approx0,
$$

则没有足够 $\beta$ reservoir。

若

$$
\Delta\nu_{\beta,i}(B)
$$

接近线性函数，则它只会 renormalize slope，不会产生强弯曲。

强 bending 要求：

$$
\frac{d}{dB}\Delta\nu_{\beta,i}(B)
$$

与

$$
g_\alpha\frac{d\phi}{dB}
$$

同量级，并且在目标场区强烈变化。

### 11.3 Compare with total LFD

最终从 total DOS 计算

$$
N(\mu,B)=\int_{-\infty}^{\mu}D(E,B)dE-N_{\rm ref}(B),
$$

反解

$$
\mu=\mu(n,B),
$$

再作

$$
D(n,B)=D(\mu(n,B),B).
$$

从 $D(n,B)$ 或 $R_{xx}$ proxy 提取 line position：

$$
\nu_i^{\rm total}(B).
$$

定义相对于 clean $\alpha$ fan 的偏移：

$$
\delta\nu_i(B)
=\nu_i^{\rm total}(B)-[\nu_C+g_\alpha i\phi(B)].
$$

若机制正确，应有

$$
\delta\nu_i(B)\sim \text{function controlled by }N_\beta.
$$

---

## 12. 失败判据

应拒绝或重新调参的情形：

1. zero-field particle side 仍是三组 disconnected $C_3$ pockets，而没有 single $C_3$-invariant $\alpha$ orbit 或 $\beta$ network；
2. $D(E,B)$ 中只有 clean LL ridges，没有 dense background；
3. $D_\beta(E,B)$ 也是清晰 LL，不是 reservoir；
4. $\Delta\nu_{\beta,i}(B)$ 很小或近似线性；
5. hole side 也出现同等强度的 $\beta$ reservoir，导致 hole fan 同样弯曲；
6. apparent bending 只来自 mask choice，而 total $D(n,B)$ 不弯；
7. mask 结果对 $\sigma_k,\sigma_E$ 极端敏感；
8. 为了得到 particle-side smaller spacing，不得不显式破坏 $C_3$ 或加入 zero-field translation-breaking inter-patch coupling。

---

## 13. 推荐实际工作流

### Step 1: 构造 zero-field $H_+(\mathbf k)$

使用

$$
H_+(\mathbf k)=d_0\sigma_0+d_z\sigma_z+d_x\sigma_x
$$

并保证 $C_3$、$C_2$ line symmetry、valley $U(1)$、TR copy。

### Step 2: 检查 Dirac cones

沿 $C_2$ line 做 continuity-based tracking 和 $C_2$ label 检查。确认 Dirac roots 和 $C_3$ copies。

### Step 3: 调 particle-hole asymmetry

通过 $d_0$ higher harmonics 调整：

$$
E<E_D:
3\text{ disconnected pockets per valley},
$$

$$
E>E_D:
1\text{ clean }C_3\text{-invariant }\alpha\text{ orbit}+\beta\text{ network}.
$$

### Step 4: 计算 zero-field DOS 和 masks

构造 $\chi_\alpha,\chi_\beta$，确认 $D_{0,\beta}$ 在 particle-side target window 中足够大。

### Step 5: 运行已有 Peierls--Hofstadter code

得到 full spectrum 和 total $D(E,B)$。

### Step 6: 计算 mask-projected DOS

使用 $P_\alpha^{(0)},P_\beta^{(0)}$ 计算

$$
D_\alpha(E,B),\quad D_\beta(E,B),
$$

以及 cumulative

$$
N_\beta(\mu,B).
$$

### Step 7: 检查 bending capacity

提取 $E_i^\alpha(B)$，计算

$$
\Delta\nu_{\beta,i}(B)
=N_\beta(E_{i+1}^\alpha,B)-N_\beta(E_i^\alpha,B).
$$

只有当该量足够大且随 $B$ 强非线性变化时，才进入最终 LFD 计算。

### Step 8: 计算 total LFD

使用 total DOS：

$$
D(E,B)\to n(\mu,B)\to \mu(n,B)\to D(n,B)\quad\text{or}\quad R_{xx}(n,B).
$$

mask-projected DOS 只用于解释，不用于替代 total DOS。

---

## 14. 参考文献与方法锚点

1. L. Onsager, “Interpretation of the de Haas-van Alphen effect,” *Philosophical Magazine* **43**, 1006--1008 (1952).
2. I. M. Lifshitz and A. M. Kosevich, “Theory of magnetic susceptibility in metals at low temperatures,” *Soviet Physics JETP* **2**, 636--645 (1956).
3. L. M. Roth, “Semiclassical theory of magnetic energy levels and magnetic susceptibility of Bloch electrons,” *Physical Review* **145**, 434--448 (1966).
4. G. H. Wannier, “A result not dependent on rationality for Bloch electrons in a magnetic field,” *Physica Status Solidi B* **88**, 757--765 (1978).
5. P. Středa, “Theory of quantised Hall conductivity in two dimensions,” *Journal of Physics C: Solid State Physics* **15**, L717--L721 (1982).
6. Pilkyung Moon, Youngwook Kim, Mikito Koshino, Takashi Taniguchi, Kenji Watanabe, and Jurgen H. Smet, “Nonlinear Landau Fan Diagram for Graphene Electrons Exposed to a Moiré Potential,” *Nano Letters* **24**, 3339--3346 (2024).

