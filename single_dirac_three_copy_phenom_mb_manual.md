# 从 Single Dirac Cone 到三份 $C_3$ Copy 的唯象 Magnetic Breakdown 模型手册

## 0. 目标与适用范围

本文档给出一个从 single Dirac cone 出发的 **phenomenological orbit-space model**。它用于在全局两能带连续模型调参困难时，快速构造和测试：

1. 三个 $C_3$-related Dirac pockets；
2. particle-hole asymmetric magnetic breakdown；
3. particle side 的 LL splitting / degeneracy lifting；
4. particle side 的 $eta$-background；
5. 因 $eta$-background 导致的 nonlinear Landau fan。

这不是 microscopic global Hamiltonian 的 exact minimal-coupling calculation。这里的 magnetic breakdown coupling 是在 semiclassical LL / orbit basis 中唯象加入的。

因此本文区分两层：

$$
\text{microscopic exact model: } H_B=H_0(\boldsymbol\Pi),
$$

和

$$
\text{phenomenological orbit model: independent LLs + MB coupling + beta background}.
$$

本文采用第二种。

---

## 1. 单个 particle-hole asymmetric Dirac cone

从一个局部 Dirac cone 出发：

$$
h_D(\mathbf q)
=
E_D\sigma_0
+
D(q)\sigma_0
+
v(q_x\sigma_x+q_y\sigma_y),
$$

其中

$$
q=|\mathbf q|.
$$

最小 particle-hole asymmetric scalar dispersion 取为

$$
D(q)=\alpha q^2.
$$

于是两支能带为

$$
E_s(q)=E_D+\alpha q^2+s v q,
\qquad
s=+1,-1.
$$

其中：

- $s=+1$ 是 particle side；
- $s=-1$ 是 hole side；
- $v$ 是 Dirac velocity；
- $\alpha$ 控制 particle-hole asymmetry。

若取

$$
\alpha<0,
$$

particle branch

$$
E_+(q)=E_D+\alpha q^2+vq
$$

向下弯曲，并在

$$
q_{\rm top}=\frac{v}{2|\alpha|}
$$

处有 band-top-like scale：

$$
E_{\rm top}=E_D+\frac{v^2}{4|\alpha|}.
$$

这会使 particle side 的 constant-energy orbits 在目标能量附近半径变大，cyclotron mass 增大，更容易发生 magnetic breakdown 和形成 dense background。

---

## 2. 三个 $C_3$-related copies

令三个 Dirac cone centers 为

$$
\mathbf K_j
=
K_0
\begin{pmatrix}
\cos(2\pi j/3)\\
\sin(2\pi j/3)
\end{pmatrix},
\qquad
j=0,1,2.
$$

相邻 centers 的距离为

$$
L=|\mathbf K_i-\mathbf K_j|=\sqrt3K_0.
$$

每个 copy 的局域动量为

$$
\mathbf q_j=R_{-2\pi j/3}(\mathbf k-\mathbf K_j).
$$

零场 copy model 是

$$
H_{\rm copy}^{(0)}(\mathbf k)
=
\bigoplus_{j=0}^{2}h_D(\mathbf q_j).
$$

这个 block-diagonal model 本身不会自动产生 MB。MB 将在 LL / orbit basis 中通过唯象 coupling 加入。

---

## 3. 未耦合的 semiclassical Landau quantization

对 circular orbit，面积为

$$
A(E)=\pi r^2(E).
$$

Dirac-like orbit 的 Onsager quantization 近似为

$$
A(E)=2\pi\frac{eB}{\hbar}N,
\qquad
N=0,1,2,\dots .
$$

因此

$$
r_N(B)=\sqrt{\frac{2eB}{\hbar}N}.
$$

未耦合 Landau levels 为

$$
E_{s,N}^{(0)}(B)
=
E_D+\alpha r_N^2+s v r_N.
$$

如果使用无量纲单位 $\hbar=e=1$，则

$$
r_N(B)=\sqrt{2BN},
$$

$$
E_{s,N}^{(0)}(B)
=
E_D+2\alpha BN+s v\sqrt{2BN}.
$$

三份 copies 在未耦合时简并：

$$
E_{s,N,j}^{(0)}(B)=E_{s,N}^{(0)}(B),
\qquad
j=0,1,2.
$$

Dirac zeroth LL $N=0$ 可以单独处理。如果目标是 particle-side nonzero LL splitting 和 fan bending，可先不对 $N=0$ 加 MB coupling，避免 central LL 干扰。

---

## 4. 从 orbit distance 定义 MB probability

对能量 $E$，设某一侧的 orbit 半径为 $r_s(E)$。在 circular approximation 下，两相邻 pockets 的最近动量距离为

$$
\Delta k_s(E)=\max\{L-2r_s(E),0\}.
$$

这里

$$
L=\sqrt3K_0.
$$

磁长为

$$
\ell_B=\sqrt{\frac{\hbar}{eB}}.
$$

MB 的快速几何判据是

$$
\Delta k_s(E)\ell_B\lesssim O(1).
$$

定义唯象 MB probability：

$$
P_s(E,B)
=
\exp[-C_s\Delta k_s^2(E)\ell_B^2].
$$

等价地，

$$
P_s(E,B)
=
\exp\left[-\frac{B_{0,s}(E)}{B}\right],
$$

其中

$$
B_{0,s}(E)=C_s\frac{\hbar}{e}\Delta k_s^2(E).
$$

$C_s$ 是 geometry factor。它吸收了 contour curvature、local velocity angle、forbidden action shape 等没有显式建模的因素。

---

## 5. Particle-hole asymmetric MB 的目标判据

目标是

$$
\text{particle side: }P_+(E_+,B_\star)\sim 0.1\text{--}1,
$$

但

$$
\text{hole side: }P_-(E_-,B_{\max})\ll1.
$$

等价地：

$$
\Delta k_+(E_+)\ell_B(B_\star)\lesssim 1,
$$

$$
\Delta k_-(E_-)\ell_B(B_{\max})\gtrsim 3.
$$

更一般地，可以写成

$$
\Delta k_+(E_+)<\frac{c_{\rm strong}}{\ell_B(B_\star)},
\qquad
c_{\rm strong}\sim1,
$$

$$
\Delta k_-(E_-)>\frac{c_{\rm weak}}{\ell_B(B_{\max})},
\qquad
c_{\rm weak}\sim3.
$$

调参上，最直接的旋钮是：

1. 取 $\alpha<0$，让 particle orbit 半径变大；
2. 调 $K_0$ 控制三个 pockets 的几何距离；
3. 调 $C_+$ 和 $C_-$ 控制 phenomenological MB threshold；
4. 让目标 particle energy $E_+$ 靠近 $E_{\rm top}$，增强 particle-side large mass / dense-level 倾向。

---

## 6. 能量到半径的反解

给定能量

$$
e=E-E_D,
$$

半径由

$$
\alpha r^2+s v r-e=0
$$

给出。

当 $\alpha\neq0$ 时，候选根为

$$
r_{s,\pm}(E)
=
\frac{-sv\pm\sqrt{v^2+4\alpha e}}{2\alpha}.
$$

物理上保留实数且非负的根。

当 $\alpha<0$ 且 $s=+$ 时，particle branch 在 $E<E_{\rm top}$ 可能有两个正根：

$$
r_{+,\rm inner}(E),
\qquad
r_{+,\rm outer}(E).
$$

一般：

- inner branch 对应 Dirac-pocket-like $\alpha$ orbit；
- outer branch 或接近 band top 的区域可作为 dense $\beta$-sector 的来源。

若只做最小模型，可先使用 inner branch 定义 $\alpha$-LL，并把 $\beta$-sector 作为额外背景 DOS 加入。

---

## 7. 在三 copy LL space 中加入 MB coupling

对每个 $(s,N)$，copy space 为

$$
|j;s,N\rangle,
\qquad
j=0,1,2.
$$

定义 MB amplitude：

$$
t_{s,N}(B)=t_{\rm scale}(E_{s,N},B)\sqrt{P_s(E_{s,N}^{(0)},B)}.
$$

注意：

$$
P_s=\text{probability},
\qquad
\sqrt{P_s}=\text{amplitude factor}.
$$

最简单取

$$
t_{\rm scale}(E,B)=t_{\rm scale}^{(0)}.
$$

更合理的选择是

$$
t_{\rm scale}(E,B)=c_t\Delta E_{\rm LL}(E,B),
$$

其中 $\Delta E_{\rm LL}$ 是附近未耦合 LL spacing，$c_t\sim0.1\text{--}1$。

构造三 copy ring Hamiltonian：

$$
H_{s,N}^{\rm copy}(B)
=
E_{s,N}^{(0)}(B)I_3
+
t_{s,N}(B)
\sum_{j=0}^{2}
\left(
 e^{i\phi_{\rm MB}/3}|j+1\rangle\langle j|
+
e^{-i\phi_{\rm MB}/3}|j\rangle\langle j+1|
\right),
$$

其中 $j+1$ 按 mod $3$ 计算。

该矩阵的 $C_3$ sectors 为

$$
m=0,1,2.
$$

本征能量为

$$
E_{s,N,m}(B)
=
E_{s,N}^{(0)}(B)
+
2t_{s,N}(B)
\cos\left(\frac{\phi_{\rm MB}+2\pi m}{3}\right).
$$

若 $\phi_{\rm MB}=0$，则得到 one singlet + one doublet：

$$
E=E^{(0)}+2t,
$$

$$
E=E^{(0)}-t
\quad
\text{twofold}.
$$

若 $\phi_{\rm MB}$ 非零且不是特殊值，则三个 $C_3$ sectors 可以完全分裂。因为外磁场破坏 time reversal，把 $\phi_{\rm MB}$ 作为 Stokes / orbit phase 是合理的 phenomenological choice。

---

## 8. 只让 particle side MB 的参数实现

定义

$$
P_+(E,B)=\exp[-B_{0,+}(E)/B],
$$

$$
P_-(E,B)=\exp[-B_{0,-}(E)/B].
$$

要求

$$
B_{0,+}(E_+)\sim B_\star,
$$

$$
B_{0,-}(E_-)\gg B_{\max}.
$$

实际实现可以写成 amplitude scale：

$$
t_{+,N}(B)
=
t_{\rm scale}
\exp\left[-\frac{B_{T,+}(E_{+,N})}{B}\right],
$$

$$
t_{-,N}(B)
=
t_{\rm scale}
\exp\left[-\frac{B_{T,-}(E_{-,N})}{B}\right],
$$

其中

$$
B_T=\frac{B_0}{2}
$$

是 amplitude exponent scale。

一个方便的 particle-side resonance choice 是

$$
B_{T,+}(E)=B_{+,\min}+\beta_+(E-E_\beta)^2,
$$

$$
B_{T,-}(E)=B_{-,0},
\qquad
B_{-,0}\gg B_{\max}.
$$

这样 particle side 在 $E\approx E_\beta$ 且 $B\sim B_\star$ 时发生明显 splitting，而 hole side 基本不 split。

---

## 9. 加入 particle-side $\beta$-background

只有三 copy LL splitting 不一定给出 fan bending。fan bending 需要 particle side 存在额外 dense states：

$$
D_\beta(E,B)\neq0.
$$

可加入 particle-side background DOS：

$$
D_\beta(E,B)
=
A_\beta P_+(E,B)D_{\beta,0}(E).
$$

一个 Lorentzian 形式为

$$
D_{\beta,0}(E)
=
\frac{1}{\pi}
\frac{\Gamma_\beta}{(E-E_\beta)^2+\Gamma_\beta^2}.
$$

也可以取 box-like dense band：

$$
D_{\beta,0}(E)
=
D_\beta^{(0)}
\Theta(E-E_{\beta,1})
\Theta(E_{\beta,2}-E).
$$

只在 particle side 打开：

$$
D_\beta(E<0,B)=0.
$$

总 DOS 写成

$$
D(E,B)=D_\alpha(E,B)+D_\beta(E,B)+D_{\rm bg}(E,B).
$$

其中 discrete $\alpha$ LL contribution 是

$$
D_\alpha(E,B)
=
\sum_{s,N,m}
d_{s,N,m}(B)
G_\Gamma(E-E_{s,N,m}(B)).
$$

如果只希望 particle-side fan bending，可以让 $D_\beta$ 只作用于 $s=+$ 的目标能区。

---

## 10. 从 spectrum / DOS 到 $R_{xx}(n,B)$

不要直接画

$$
n=n_0+\nu B.
$$

应走统一流程：

$$
E_j(B)
\rightarrow
D(E,B)
\rightarrow
n(\mu,B)
\rightarrow
\mu(n,B)
\rightarrow
R_{xx}(n,B).
$$

累计 density：

$$
n(\mu,B)=\int_{E_{\rm min}}^{\mu}D(E,B)dE-n_{\rm ref}(B).
$$

若相邻 $\alpha$-LL 中间有 $\beta$ background，则填充增量变成

$$
\Delta n
=
g_\alpha\frac{eB}{h}
+
\int_{E_i^\alpha}^{E_{i+1}^\alpha}D_\beta(E,B)dE.
$$

当 $D_\beta=0$ 时，fan line 近似线性；当 $D_\beta\neq0$ 时，fan line 可以弯曲。

---

## 11. 参数调节流程

### Step 1：固定 single cone 参数

先取无量纲单位：

$$
v=1,
\qquad
E_D=0.
$$

取

$$
\alpha<0
$$

制造 particle-side flattening。

推荐初始：

$$
\alpha=-0.15\text{--}-0.35.
$$

对应 band-top scale：

$$
E_{\rm top}=\frac{1}{4|\alpha|}.
$$

若希望 particle-side bending 出现在 $E_\beta$，应让

$$
E_\beta\lesssim E_{\rm top}.
$$

### Step 2：调 $K_0$ 设定 MB 磁场尺度

中心距离为

$$
L=\sqrt3K_0.
$$

$K_0$ 越大，pockets 越远，MB 越晚发生；$K_0$ 越小，MB 越早发生。

用

$$
\Delta k_+(E_+)\ell_B(B_\star)\sim1
$$

确定初始 $K_0$。

### Step 3：调 $C_+$ 与 $C_-$

如果希望只做几何驱动，可以取

$$
C_+=C_-=1.
$$

如果需要更强 phenomenological asymmetry，可以取

$$
C_+\sim0.5\text{--}1.5,
$$

$$
C_-\sim5\text{--}20.
$$

这等价于让 hole-side forbidden action 更大。

### Step 4：调 splitting strength

选择

$$
t_{\rm scale}=0.02\text{--}0.10.
$$

要求 particle-side splitting 大于 DOS broadening：

$$
2|t_{+,N}(B)|\gtrsim \Gamma.
$$

hole side 要求

$$
2|t_{-,N}(B)|\ll \Gamma.
$$

### Step 5：加入 $\beta$ background

选择

$$
E_\beta\sim0.3\text{--}0.6,
$$

$$
\Gamma_\beta\sim0.05\text{--}0.2,
$$

$$
A_\beta\sim0.5\text{--}3.
$$

若 fan bending 不明显，增大 $A_\beta$ 或增大 $\Gamma_\beta$；若 particle side LL 被完全洗掉，减小 $A_\beta$ 或 $\Gamma_\beta$。

---

## 12. 检验标准

### 12.1 Particle-side MB

要求

$$
P_+(E_{+,N},B_\star)\sim0.1\text{--}1.
$$

对应 splitting

$$
\Delta E_{\rm split}
\sim
2t_{\rm scale}\sqrt{P_+}
$$

应大于可视化展宽 $\Gamma$。

### 12.2 Hole-side no MB

要求

$$
P_-(E_{-,N},B_{\max})<10^{-2},
$$

或至少

$$
2t_{\rm scale}\sqrt{P_-}\ll\Gamma.
$$

hole side 应保持近似 unresolved three-copy degeneracy 或普通 clean LL ladder。

### 12.3 Beta background

要求 particle side 有

$$
D_\beta(E,B)>0
$$

且相邻 $\alpha$-LL 之间的 integrated beta density 不可忽略：

$$
\int_{E_i^\alpha}^{E_{i+1}^\alpha}D_\beta(E,B)dE
\sim
O\left(g_\alpha\frac{eB}{h}\right)
$$

或至少达到其明显 fraction。

### 12.4 Fan bending

最终用 DOS pipeline 检查：

$$
D(E,B)
\rightarrow
n(\mu,B)
\rightarrow
D(n,B).
$$

合格表现：

1. particle side fan line 在目标磁场区间弯曲；
2. hole side fan line 保持近似线性；
3. particle side 有 dense background；
4. hole side 没有对应 dense background。

---

## 13. Python 实现骨架

```python
import numpy as np


def positive_roots_quadratic(a, b, c, tol=1e-12):
    disc = b*b - 4*a*c
    if disc < -tol:
        return []
    disc = max(disc, 0.0)
    if abs(a) < tol:
        if abs(b) < tol:
            return []
        r = -c / b
        return [r] if r >= -tol else []
    roots = [(-b + np.sqrt(disc))/(2*a), (-b - np.sqrt(disc))/(2*a)]
    return sorted([r for r in roots if r >= -tol])


def orbit_radii(E, side, ED=0.0, v=1.0, alpha=-0.2):
    e = E - ED
    s = +1 if side == "particle" else -1
    # alpha r^2 + s v r - e = 0
    return positive_roots_quadratic(alpha, s*v, -e)


def choose_alpha_radius(E, side, ED=0.0, v=1.0, alpha=-0.2):
    roots = orbit_radii(E, side, ED=ED, v=v, alpha=alpha)
    if len(roots) == 0:
        return np.nan
    # inner branch is the Dirac-like alpha orbit
    return roots[0]


def delta_k(E, side, K0=1.0, ED=0.0, v=1.0, alpha=-0.2):
    L = np.sqrt(3.0) * K0
    r = choose_alpha_radius(E, side, ED=ED, v=v, alpha=alpha)
    if not np.isfinite(r):
        return np.nan
    return max(L - 2.0*r, 0.0)


def P_MB(E, B, side, K0=1.0, ED=0.0, v=1.0, alpha=-0.2, Cgeom=1.0):
    dk = delta_k(E, side, K0=K0, ED=ED, v=v, alpha=alpha)
    if not np.isfinite(dk):
        return 0.0
    ellB = 1.0 / np.sqrt(B)  # hbar=e=1
    return np.exp(-Cgeom * (dk * ellB)**2)


def dirac_LL_energy(N, B, side, ED=0.0, v=1.0, alpha=-0.2):
    rN = np.sqrt(2.0 * B * N)  # hbar=e=1, Dirac gamma=0
    s = +1 if side == "particle" else -1
    return ED + alpha * rN**2 + s * v * rN


def split_three_copy_levels(
    B,
    Nmax=30,
    side="particle",
    ED=0.0,
    v=1.0,
    alpha=-0.2,
    K0=1.0,
    t_scale=0.05,
    Cgeom=1.0,
    phi_MB=0.4,
):
    levels = []
    meta = []

    for N in range(1, Nmax + 1):
        E0 = dirac_LL_energy(N, B, side=side, ED=ED, v=v, alpha=alpha)
        P = P_MB(abs(E0-ED), B, side=side, K0=K0, ED=0.0, v=v, alpha=alpha, Cgeom=Cgeom)
        t = t_scale * np.sqrt(P)

        for m in range(3):
            Em = E0 + 2.0 * t * np.cos((phi_MB + 2.0*np.pi*m)/3.0)
            levels.append(Em)
            meta.append((side, N, m, P, t))

    return np.array(levels), meta


def gaussian(Egrid, E0, Gamma):
    return np.exp(-0.5*((Egrid-E0)/Gamma)**2)/(np.sqrt(2*np.pi)*Gamma)


def beta_dos(Egrid, B, Ebeta=0.4, Gbeta=0.08, Abeta=1.0, Pbeta=1.0):
    lor = (1.0/np.pi) * Gbeta / ((Egrid - Ebeta)**2 + Gbeta**2)
    return Abeta * Pbeta * lor


def total_dos(Egrid, B, params):
    Gamma = params.get("Gamma_LL", 0.03)
    Nmax = params.get("Nmax", 30)
    t_scale = params.get("t_scale", 0.05)
    K0 = params.get("K0", 1.0)
    v = params.get("v", 1.0)
    alpha = params.get("alpha", -0.2)
    ED = params.get("ED", 0.0)
    phi_MB = params.get("phi_MB", 0.4)
    Cp = params.get("C_particle", 1.0)
    Ch = params.get("C_hole", 10.0)

    dos = np.zeros_like(Egrid)

    for side, Cgeom in [("particle", Cp), ("hole", Ch)]:
        levels, meta = split_three_copy_levels(
            B,
            Nmax=Nmax,
            side=side,
            ED=ED,
            v=v,
            alpha=alpha,
            K0=K0,
            t_scale=t_scale,
            Cgeom=Cgeom,
            phi_MB=phi_MB,
        )
        for E in levels:
            dos += B * gaussian(Egrid, E, Gamma)

    # particle-side beta background
    Ebeta = params.get("Ebeta", 0.4)
    Gbeta = params.get("Gbeta", 0.08)
    Abeta = params.get("Abeta", 1.0)

    Pbeta = P_MB(Ebeta, B, side="particle", K0=K0, ED=ED, v=v, alpha=alpha, Cgeom=Cp)
    dos += beta_dos(Egrid, B, Ebeta=Ebeta, Gbeta=Gbeta, Abeta=Abeta, Pbeta=Pbeta)

    return dos
```

---

## 14. 推荐初始参数

一组可直接测试的无量纲参数是：

```python
params = dict(
    ED=0.0,
    v=1.0,
    alpha=-0.25,
    K0=1.0,
    Nmax=40,
    Gamma_LL=0.03,
    t_scale=0.06,
    phi_MB=0.5,
    C_particle=1.0,
    C_hole=12.0,
    Ebeta=0.45,
    Gbeta=0.10,
    Abeta=1.5,
)
```

需要扫描的主参数是：

$$
\alpha,
\quad
K_0,
\quad
C_+,
\quad
C_-,
\quad
t_{\rm scale},
\quad
A_\beta,
\quad
E_\beta,
\quad
\Gamma_\beta.
$$

---

## 15. 常见失败模式

### 15.1 particle side 不发生 MB

表现：

$$
P_+(E,B_\star)\ll1.
$$

修正：

1. 减小 $K_0$；
2. 增大 $|\alpha|$；
3. 减小 $C_+$；
4. 增大目标 $E_\beta$ 靠近 band top；
5. 增大 $t_{\rm scale}$。

### 15.2 hole side 也发生 MB

表现：

$$
P_-(E,B_{\max})\not\ll1.
$$

修正：

1. 增大 $C_-$；
2. 增大 $K_0$；
3. 减小 $|\alpha|$；
4. 限制 $D_\beta(E<0,B)=0$；
5. 减小 hole-side $t_-$ 或直接设 $t_-=0$ 作为 phenomenological limit。

### 15.3 fan bending 不明显

表现：particle LL splitting 有了，但 $D(n,B)$ 仍近似线性。

修正：

1. 增大 $A_\beta$；
2. 增大 $\Gamma_\beta$；
3. 把 $E_\beta$ 放在相邻 $\alpha$-LL 中间；
4. 使 $D_\beta$ 随 $B$ 增强，例如乘 $P_+(E,B)$；
5. 检查是否真正通过 $D(E,B)\to n(\mu,B)$ 反演，而不是直接画 LL lines。

### 15.4 discrete LL 被完全洗掉

表现：particle side 没有清晰 $\alpha$ ridges。

修正：

1. 减小 $A_\beta$；
2. 减小 $\Gamma_\beta$；
3. 减小 $t_{\rm scale}$；
4. 保留一部分 isolated $\alpha$-LL weight。

---

## 16. 最终合格标准

一个可用参数组应满足：

1. particle side:

$$
P_+(E_+,B_\star)\sim0.1\text{--}1;
$$

2. hole side:

$$
P_-(E_-,B_{\max})<10^{-2};
$$

3. particle-side splitting:

$$
2t_{\rm scale}\sqrt{P_+}\gtrsim \Gamma_{\rm LL};
$$

4. hole-side splitting:

$$
2t_{\rm scale}\sqrt{P_-}\ll \Gamma_{\rm LL};
$$

5. beta background:

$$
\int_{E_i^\alpha}^{E_{i+1}^\alpha}D_\beta(E,B)dE
$$

不可忽略；

6. final DOS/filling inversion 中：

$$
\text{particle fan bends, hole fan remains nearly linear}.
$$

---

## 17. 结论

当全局两能带 continuous Hamiltonian 太难调参时，可以先使用这个 three-copy phenomenological MB model。

它的物理结构是：

$$
\text{single PH-asymmetric Dirac cone}
\rightarrow
3\ C_3\text{ copies}
\rightarrow
\text{phenomenological MB splitting}
\rightarrow
\text{particle-side }\beta\text{ background}
\rightarrow
D(E,B)\rightarrow R_{xx}(n,B).
$$

这套方法牺牲 microscopic exactness，但参数直接、稳定、容易扫描。它适合作为寻找目标 Landau fan 结构的第一阶段模型，然后再把有效机制反推回 global continuous Hamiltonian。
