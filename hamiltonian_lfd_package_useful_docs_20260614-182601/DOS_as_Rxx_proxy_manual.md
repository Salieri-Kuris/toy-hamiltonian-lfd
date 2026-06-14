# DOS 作为 $R_{xx}$ proxy 的唯象手册

## 0. 目的与适用范围

本手册用于如下问题：已经通过 Hofstadter / Peierls 磁超胞严格求解得到总 DOS

$$
D_{\rm Hof}(E,B),
$$

但暂时不计算完整 Kubo conductivity，希望从 DOS 唯象地构造一个可用于比较 Landau fan diagram 的

$$
R_{xx}^{\rm proxy}(n,B).
$$

核心目标不是定量复现实验电阻幅值，而是定性判断：

1. 哪些 Hofstadter gaps / low-DOS regions 可能对应实验中的 $R_{xx}$ minima；
2. 哪些 DOS ridges 可能对应 $R_{xx}$ maxima 或 $\sigma_{xx}$ maxima；
3. particle-hole asymmetry、background DOS、finite broadening、density inversion 如何改变 fan-line positions 和 visibility。

本手册不引入额外的 magnetic-breakdown visibility factor。若 Hofstadter spectrum 已经由完整磁场哈密顿量给出，则 magnetic-field-induced hybridization 已经包含在谱内；后处理只允许加入 disorder、temperature、density inhomogeneity、mobility / extended-state filtering 等输运可见性参数。

---

## 1. 基本原则

### 1.1 总 DOS 决定 density

总 density 必须由总 DOS 累计得到：

$$
n(\mu,B)
=
\int_{-\infty}^{\mu}D_{\rm tot}(E,B)\,dE
-
n_{\rm ref}(B).
$$

这里 $D_{\rm tot}$ 应包括所有 states。不能为了隐藏某个 pocket 的 fan line 而把它从 $D_{\rm tot}$ 中删掉。否则它也不会贡献累计态数，也就不能改变 $\mu(n,B)$ 或扭曲其他 fan lines。

### 1.2 $R_{xx}$ proxy 更接近 visible / extended DOS

真实 $R_{xx}$ 不是总 DOS 的函数。它还取决于 extended states、mobility edge、scattering time、Hall tensor、edge states、contacts 和 interaction-induced gaps。唯象层面可以写为

$$
R_{xx}^{\rm proxy}(n,B)
\propto
D_{\rm ext}(\mu(n,B),B),
$$

其中 $D_{\rm ext}$ 是 extended-state DOS 或 visible DOS 的近似。

因此最重要的区分是：

$$
D_{\rm tot} \quad \text{用于 } n(\mu,B),
$$

$$
D_{\rm vis}\text{ 或 }D_{\rm ext} \quad \text{用于 } R_{xx}^{\rm proxy}.
$$

### 1.3 density inversion 是必要步骤

实验横轴通常是 density / filling，而不是 chemical potential。因此必须计算

$$
D(E,B)
\longrightarrow
n(\mu,B)
\longrightarrow
\mu(n,B)
\longrightarrow
R_{xx}^{\rm proxy}(n,B).
$$

不要先在 $n$-$B$ 平面中假设 fan lines 再画 Gaussian ridge / dip。

---

## 2. 输入数据

假设 Hofstadter 求解给出离散磁场 $B_i$ 上的能谱或 DOS。

### 2.1 若输入是能级 / 磁子带中心

若有能级

$$
\epsilon_\lambda(B),
$$

以及权重

$$
w_\lambda(B),
$$

则可以直接构造 broadened DOS。

### 2.2 若输入已经是网格 DOS

若已经有

$$
D_{\rm Hof}(E_a,B_i),
$$

则可以在能量方向做归一化卷积，得到 broadened DOS。卷积核必须保持总态数守恒：

$$
\int dE\,D_\Gamma(E,B)
=
\int dE\,D_{\rm Hof}(E,B).
$$

---

# 方法 A+B：直接 broadened-DOS proxy

方法 A+B 是最小且最安全的方案。它不把总 DOS 人为分成 $\alpha,\beta$ sectors，而是用全局或标量能量依赖的 broadening 处理整个 DOS。

---

## 3. 方法 A：全局 / 磁场 / 温度 broadening

### 3.1 定义

对每个 Hofstadter level 使用 Gaussian broadening：

$$
D_\Gamma(E,B)
=
\sum_\lambda
w_\lambda(B)
G_{\Gamma_\lambda(B)}(E-\epsilon_\lambda(B)),
$$

其中

$$
G_\Gamma(x)
=
\frac{1}{\sqrt{2\pi}\Gamma}
\exp\left[-\frac{x^2}{2\Gamma^2}\right].
$$

最小模型取

$$
\Gamma_\lambda(B)=\Gamma_0.
$$

更一般地取

$$
\Gamma_{\rm eff}^2(B)
=
\Gamma_0^2
+
\Gamma_T^2
+
\Gamma_B^2(B).
$$

常用磁场依赖为

$$
\Gamma_B(B)=a\sqrt{B},
$$

或

$$
\Gamma_B(B)=bB.
$$

有限温度可以先用

$$
\Gamma_T=c_T k_B T
$$

近似并入 effective broadening。严格有限温处理应使用 Fermi function 卷积，但对于定性 LFD 可视化，$\Gamma_T$ 近似通常足够。

### 3.2 参数物理含义

| 参数 | 物理含义 | 调节效果 |
|---|---|---|
| $\Gamma_0$ | 残余 disorder / quantum lifetime broadening / 数值分辨率 | 洗掉小 Hofstadter gaps；控制全局线宽 |
| $c_T k_BT$ | finite-temperature smearing | 温度越高，fan lines 越宽，弱 minima 越不清楚 |
| $a\sqrt B$ | Landau-level disorder broadening 的常见唯象形式 | 高场 LL 变宽；小 gaps 更容易消失 |
| $bB$ | 更强的 field-dependent broadening | 高场细结构迅速被洗掉 |

### 3.3 参数如何选取

第一步只使用全局 $\Gamma_0$。令 $\Delta_{\rm vis}$ 是希望实验可见的主要 gap，$\Delta_{\rm inv}$ 是希望被洗掉的小 gap。合理窗口是

$$
\Delta_{\rm inv}
\lesssim
\Gamma_0
\lesssim
0.2\text{--}0.4\,\Delta_{\rm vis}.
$$

如果 $\Gamma_0$ 大到会洗掉 hole-side 清楚的 sixfold fan，则模型失败，不能靠 smearing 解释目标图像。

若实验高场线宽明显变大，再加入 $a\sqrt B$ 或 $bB$。不要一开始就同时调多个 broadening 参数。

### 3.4 从 broadened DOS 到 density

在每个固定 $B$ 上计算

$$
n(\mu,B)
=
\int_{E_{\min}}^{\mu}D_\Gamma(E,B)\,dE
-
n(\mu_{\rm ref},B).
$$

数值上用 cumulative trapezoidal integral 即可。

然后反解

$$
\mu=\mu(n,B).
$$

最后得到

$$
D_\Gamma(n,B)
=
D_\Gamma(\mu(n,B),B).
$$



---

## 4. 方法 B：标量能量依赖 broadening

方法 B 是方法 A 的物理增强版。它允许不同能量区域有不同 smearing，但要求 broadening 来自全局标量规则，而不是手工指定“这个区域是 $\alpha$、那个区域是 $\beta$”。

### 4.1 推荐形式

取

$$
\Gamma(E,B)
=
\sqrt{
\Gamma_0^2
+
\Gamma_T^2
+
\Gamma_B^2(B)
+
\Gamma_D^2(E)
}.
$$

最有物理依据的标量选择是

$$
\Gamma_D(E)
=
\lambda_D\,\overline D_0(E),
$$

其中 $\overline D_0(E)$ 是归一化后的零场 DOS。也可以在已知 vHS 能量时加入

$$
\Gamma_{\rm vH}(E)
=
\Gamma_{\rm vH}^{0}
\exp\left[-\frac{(E-E_{\rm vH})^2}{2w_{\rm vH}^2}\right].
$$

此时

$$
\Gamma^2(E,B)
=
\Gamma_0^2+
\Gamma_T^2+
\Gamma_B^2(B)
+
\lambda_D^2\overline D_0^2(E)
+
\Gamma_{\rm vH}^2(E).
$$

### 4.2 为什么合理

在 Born approximation 中，散射率通常与可散射终态数有关：

$$
\frac{1}{\tau_q(E)}\propto D_0(E),
$$

因此

$$
\Gamma(E)=\frac{\hbar}{2\tau_q(E)}
$$

可以随 $D_0(E)$ 增大而增大。

物理后果是：

- Dirac-like 低 DOS 区域：$\Gamma(E)$ 小，LL / Hofstadter gaps 更清楚；
- vHS、flat miniband、heavy pocket 区域：$\Gamma(E)$ 大，细小 oscillations 被洗掉；
- particle-hole asymmetric band structure 会自然导致 particle side 与 hole side 的 visibility 不同。

### 4.3 如何保持态数守恒

如果 $\Gamma$ 依赖能量，推荐让每个 level 使用其中心能量处的宽度：

$$
\Gamma_\lambda(B)
=
\Gamma(\epsilon_\lambda(B),B).
$$

然后写

$$
D_\Gamma(E,B)
=
\sum_\lambda
w_\lambda(B)
G_{\Gamma(\epsilon_\lambda(B),B)}(E-\epsilon_\lambda(B)).
$$

不要写成

$$
G_{\Gamma(E,B)}(E-\epsilon_\lambda),
$$

因为这通常会破坏单个 peak 的归一化，从而人为改变累计 density。

### 4.4 参数如何选取

1. 先取 $\lambda_D=0$，只用方法 A 的均匀 broadening。
2. 若 particle side 确有零场 high-DOS / vHS / heavy-pocket 区域，再打开 $\lambda_D$。
3. $\lambda_D$ 的大小应满足：

$$
\Gamma(E_h,B)\ll \Delta_h^{\rm LL},
$$

但

$$
\Gamma(E_p,B)\gtrsim \Delta_p^{\rm fine},
$$

其中 $E_h$ 是 hole-side 清晰 fan 的能量窗口，$E_p$ 是 particle-side dense pocket / vHS / heavy pocket 的能量窗口。

4. 若必须取很大的 $\lambda_D$ 才能隐藏 particle-side unwanted fan，但该 $\lambda_D$ 同时洗掉 hole-side sixfold fan，则模型失败。

### 4.5 方法 A+B 如何帮助产生目标 LFD

目标 LFD 要求：

$$
\text{hole side: clear }g=6\text{ fan},
$$

$$
\text{particle side: dense background + visible fewer-fold fan + strong bending}.
$$

方法 A+B 的作用如下。

#### Hole side

如果 hole side Dirac LL spacing 大、零场 DOS 低，则

$$
\Gamma(E_h,B)\ll \Delta_h^{\rm LL}.
$$

因此 sixfold Dirac fan 的 minima 保留。若 lower branch 由于 lattice periodicity 形成额外 pocket，但该 pocket 的 fine gaps 小于 broadening，或者其 smooth cumulative DOS 的 $B$-依赖弱，则它不会形成清晰 fan，也不会明显弯曲 lower Dirac fan。

#### Particle side

如果 particle side 靠近 vHS / heavy pocket / miniband overlap，则

$$
D_0(E_p)\text{ large}
\quad\Rightarrow\quad
\Gamma(E_p,B)\text{ large}.
$$

于是 dense pocket 的细小 oscillations 不清楚，只留下 smooth DOS background。这个 background 仍然进入

$$
n(\mu,B)=\int^{\mu}D_{\rm tot}(E,B)dE,
$$

所以它可以通过 density inversion 改变可见 upper-Dirac fan 的 $n$-位置。

可见 fan line 的位置近似满足

$$
n_m(B)
=
n_C
+
 n_m^{\rm vis}(B)
+
 n_{\rm bg}(\mu_m(B),B).
$$

其斜率为

$$
\frac{dn_m}{dB}
=
\frac{dn_m^{\rm vis}}{dB}
+
\partial_B n_{\rm bg}
+
D_{\rm bg}\frac{d\mu_m}{dB}.
$$

若背景项足够大且符号与可见 fan 的裸斜率相反，则 particle-side fan 可以强烈弯曲，甚至方向反转。

---

# 方法 C：smooth background + oscillatory visibility 分离

方法 C 适合处理如下情况：总 DOS 中已经存在 dense background 和 oscillatory ridges，但总 DOS 本身没有清晰的 state labels。此时不直接对整个 DOS 解释为 $R_{xx}$，而是把总 DOS 分为平滑背景和振荡部分。

---

## 5. 方法 C 的定义

先对总 DOS 做平滑分解：

$$
D_{\rm tot}(E,B)
=
\overline D(E,B)+\delta D(E,B).
$$

其中：

- $\overline D(E,B)$：在能量方向、磁场方向或 $1/B$ 方向平滑后的 background；
- $\delta D(E,B)$：Landau / Hofstadter oscillatory component。

然后定义 visible DOS：

$$
D_{\rm vis}(E,B)
=
\overline D(E,B)
+
A_{\rm osc}(E,B)\delta D(E,B).
$$

这里 $A_{\rm osc}$ 只控制 oscillation visibility，不改变总态数。density inversion 仍然使用

$$
D_{\rm tot}(E,B),
$$

而 $R_{xx}$ proxy 使用

$$
D_{\rm vis}(E,B).
$$

---

## 6. Oscillation visibility 参数

推荐使用 thermal damping 与 Dingle damping：

$$
A_{\rm osc}(E,B)
=
R_T(E,B)R_D(E,B).
$$

其中

$$
R_T
=
\frac{X}{\sinh X},
\qquad
X=
\frac{2\pi^2 k_BT}{\hbar\omega_c(E,B)},
$$

$$
R_D
=
\exp\left[-\frac{\pi}{\omega_c(E,B)\tau_q(E)}\right].
$$

cyclotron frequency 可写为

$$
\omega_c(E,B)=\frac{eB}{m_c(E)}.
$$

cyclotron mass 可由零场 Fermi contour 面积得到：

$$
m_c(E)
=
\frac{\hbar^2}{2\pi}\frac{\partial A(E)}{\partial E}.
$$

若没有可靠的 Fermi contour 数据，也可以从 Hofstadter DOS 中局部 LL spacing 估计：

$$
\hbar\omega_c(E,B)
\sim
\Delta E_{\rm local}(E,B).
$$

---

## 7. 方法 C 的物理含义

方法 C 的核心是：

$$
\text{smooth DOS background contributes to density},
$$

但

$$
\text{only sufficiently visible oscillatory DOS contributes clear fan lines}.
$$

这比简单 $R_{xx}\propto D_{\rm tot}$ 更接近实验。因为实验中的 $R_{xx}$ minima/ridges 通常来自 extended states 和可分辨 Landau oscillations，而 smooth background 主要改变 chemical potential、quantum capacitance 和整体电导背景。

如果某个 pocket 很重，则

$$
m_c(E)\text{ large},
\qquad
\omega_c(E,B)\text{ small},
$$

于是

$$
R_T\ll1,
\qquad
R_D\ll1.
$$

它的 oscillatory fan line 不明显，但它的 total DOS 仍然进入 $n(\mu,B)$。这正是“pocket 不形成自己的清晰 fan line，却能作为 continuous background 扭曲其他 fan line”的物理机制。

---

## 8. 方法 C 如何帮助产生目标 LFD

### 8.1 Hole side 保持清晰

若 hole-side Dirac sector 有较小 cyclotron mass 或较大的 LL spacing，则

$$
\omega_c^h\text{ large},
\qquad
R_T^hR_D^h\approx1.
$$

因此

$$
D_{\rm vis}^h\approx D_{\rm tot}^h,
$$

清晰的 sixfold fan 保留。

### 8.2 Hole-side extra pocket 不污染 fan

若 lower-branch extra pocket 的 oscillation amplitude 小，或其 LL spacing 小于 disorder / temperature scale，则

$$
R_T^{p,h}R_D^{p,h}\ll1.
$$

它不形成清晰 fan line。若其 smooth cumulative density 的 $B$-依赖也弱，则它不会显著弯曲 lower Dirac fan。

### 8.3 Particle-side dense pocket 给连续背景

若 particle-side pocket 位于 high-DOS / heavy-mass / near-vHS 区域，则

$$
R_T^{p}R_D^{p}\ll1,
$$

所以 pocket 的 oscillatory component 被压制：

$$
A_{\rm osc}\delta D_p\approx0.
$$

但它的 background $\overline D_p$ 仍然保留，并进入 density inversion：

$$
n(\mu,B)=\int^\mu D_{\rm tot}(E,B)dE.
$$

因此它可以使 upper Dirac fan 的 density position 变为

$$
n_m(B)
=
n_m^{\rm upper}(B)+n_{\rm background}(\mu_m(B),B).
$$

当 background contribution 的 $B$-dependence 足够强时，fan-line slope 可以被明显改变。



---

## 11. 推荐实施顺序

### Step 1：uniform broadening baseline

使用

$$
\Gamma(E,B)=\Gamma_0.
$$

得到

$$
D_\Gamma(E,B),
\quad
n(\mu,B),
\quad
\mu(n,B),
\quad
R_{xx}^{\rm DOS}(n,B).
$$

目的：确认谱本身是否已经含有需要的 fan-line structure。

### Step 2：加入 temperature / field broadening

使用

$$
\Gamma^2(B)=\Gamma_0^2+(c_Tk_BT)^2+(a\sqrt B)^2.
$$

目的：模拟实验分辨率和高场 broadening。

### Step 3：加入 scalar energy-dependent broadening

使用

$$
\Gamma^2(E,B)=
\Gamma_0^2+
(c_Tk_BT)^2+
(a\sqrt B)^2+
\lambda_D^2\overline D_0^2(E).
$$

目的：让 zero-field high-DOS / vHS / heavy-pocket 区域自然更宽，避免手工指定 $\alpha,\beta$ labels。

### Step 4：方法 C 独立检查

构造

$$
D_{\rm tot}=\overline D+\delta D,
$$

以及

$$
D_{\rm vis}=\overline D+R_TR_D\delta D.
$$

用 $D_{\rm tot}$ 做 density inversion，用 $D_{\rm vis}$ 做 $R_{xx}$ proxy。

目的：检验“continuous background changes density but does not create visible fan lines”的机制是否成立。

### Step 5：输出两类图

1. $D_{\rm tot}(\mu(n,B),B)$：看总谱 texture；
2. $R_{xx}^{\rm proxy}(n,B)$：看实验可见 fan minima / ridges。

若实验横轴是 gate voltage，还输出 $R_{xx}^{\rm proxy}(V_g,B)$。

---

## 12. 成功与失败标准

### 成功标准

一个合理的 DOS-proxy 模型应满足：

1. uniform broadening 下已经存在主要 spectral tendency；
2. hole side 的 $g=6$ fan 在合理 $\Gamma$ 下保持清晰；
3. particle side high-DOS / heavy-pocket region 在 scalar broadening 或方法 C 下变成 smooth background；
4. smooth background 通过 density inversion 改变 upper-side fan positions；
5. 所有参数由全局物理规则定义，而不是按目标 fan line 局部调参。

### 失败标准

应拒绝模型，如果出现以下情况：

1. 只有人为指定某个 $n$-$B$ 区域的大 smearing 才能得到目标弯曲；
2. 所需 $\Gamma$ 会同时洗掉实验中清晰的 hole-side sixfold fan；
3. total DOS 中没有足够 background，但后处理仍试图强行扭曲 fan line；
4. density inversion 使用的 DOS 与 $R_{xx}$ proxy 使用的 DOS 混淆，导致态数不守恒；
5. 不同能量区域的 smearing 没有来自 $D_0(E)$、$m_c(E)$、temperature、disorder 或 density inhomogeneity 等物理输入。

---

## 13. 最简物理总结

方法 A+B 的物理图像是：

$$
\text{Hofstadter DOS}
+
\text{disorder / temperature broadening}
+
\text{density inversion}
\Rightarrow
R_{xx}^{\rm proxy}.
$$

它适合建立 baseline，并用 zero-field DOS 或 vHS 信息引入物理可解释的 energy-dependent visibility。

方法 C 的物理图像是：

$$
\text{total DOS controls density},
$$

但

$$
\text{only visible oscillatory DOS controls clear }R_{xx}\text{ fan lines}.
$$

因此它最适合表达：

$$
\text{pocket contributes continuous background}
$$

但

$$
\text{pocket does not form its own visible fan line}.
$$

对于目标 LFD，合理参数区间应实现：

$$
\text{hole side: low DOS, large LL spacing, weak damping}
\Rightarrow
\text{clear }6\text{-fold fan},
$$

$$
\text{particle side: high DOS / heavy pocket / vHS, strong damping}
\Rightarrow
\text{continuous background},
$$

$$
\text{visible upper Dirac fan}+
\text{strong background density inversion}
\Rightarrow
\text{bent or reversed particle-side fan lines}.
$$

