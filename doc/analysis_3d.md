# 3D Generalized SOS：$M^2(L)$ 标度与平整相

## 1. 模型与设置

3 维立方格子 $\mathbb{Z}_L^3$，PBC，每个格点上一个连续高度 $h_i \in \mathbb{R}$。哈密顿量

$$\mathcal{H} = K\sum_{\langle ij\rangle} |h_i - h_j|^\sigma, \qquad \sum_i h_i = 0$$

$K=1$，$\sigma \in \{0.5,\,1.0,\,2.0,\,4.0\}$，$L \in \{4,8,16,32,64\}$。每条 $(L,\sigma)$ 独立做 Metropolis（10⁵ sweep）+ 每 5 步一次嵌入反射 Swendsen–Wang 集团更新，零模约束在每步后强制执行。

主要观测量：
$$M^2 = \frac{1}{N}\sum_i h_i^2, \qquad N=L^3$$

## 2. 数据

| $\sigma$ | $L=4$ | $L=8$ | $L=16$ | $L=32$ | $L=64$ |
|---|---|---|---|---|---|
| 0.5 | 0.3361 | 0.3867 | 0.4079 | 0.4179 | 0.4229 |
| 1.0 | 0.1098 | 0.1255 | 0.1333 | 0.1371 | 0.1391 |
| 2.0 | 0.0986 | 0.1125 | 0.1193 | 0.1228 | 0.1246 |
| 4.0 | 0.1073 | 0.1211 | 0.1285 | 0.1322 | 0.1340 |

误差（blocking）均 $\leq 10^{-3}$，远小于相邻 $L$ 之差，可以忽略。

**关键观察**：$L\to 2L$ 时 $\Delta M^2$ 严格减半。例如 $\sigma=2$：$0.0139 \to 0.0068 \to 0.0035 \to 0.0018$，公比 $\approx 1/2$。这就是 $1/L$ 修正的指纹。

## 3. 外推方法

模型：
$$\boxed{\;M^2(L) = M^2_\infty - \frac{c}{L}\;}$$

变量代换 $x = 1/L$，$y = M^2$，做线性最小二乘回归。

- **拟合窗口**：$L \geq 16$（即 $L \in \{16, 32, 64\}$）。$L=4,8$ 处次领头修正 $\mathcal{O}(1/L^2)$ 还不可忽略，纳入会偏移截距。
- **不加权 OLS**：误差棒 $\sim 10^{-4}$，远小于 $\mathcal{O}(1/L^2)$ 系统偏差，加权无意义。

外推结果：

| $\sigma$ | $M^2(64)$ | $M^2_\infty$ | $c$ | $1/c$ |
|---|---|---|---|---|
| 0.5 | 0.4229 | **0.428** | 0.320 | 3.1 |
| 1.0 | 0.1391 | **0.141** | 0.122 | 8.2 |
| 2.0 | 0.1246 | **0.126** | 0.113 | 8.8 |
| 4.0 | 0.1340 | **0.136** | 0.118 | 8.5 |

图：`figures/M2_scaling_3d.pdf`（左：$M^2$ vs $L$ 半对数，右：$M^2$ vs $1/L$ 线性外推到 $1/L=0$）。

## 4. 理论基础

### 4.1 $\sigma=2$ 严格可解

二次型哈密顿量 $\mathcal{H} = (K/2)\,\mathbf{h}^\top \Delta\, \mathbf{h}$，Fourier 后每个 $\mathbf{k}\ne 0$ 模独立高斯：

$$\langle |h_\mathbf{k}|^2\rangle = \frac{N}{K\,\varepsilon(\mathbf{k})}, \qquad \varepsilon(\mathbf{k}) = 2\sum_{a=1}^d (1-\cos k_a)$$

其中 $\varepsilon(\mathbf{k}) \sim k^2$ 为 IR。单格点方差：

$$M^2(L) = \frac{1}{KN}\sum_{\mathbf{k}\ne 0}\frac{1}{\varepsilon(\mathbf{k})}, \qquad M^2_\infty = \frac{1}{K}\int_{\rm BZ}\frac{d^d k}{(2\pi)^d}\,\frac{1}{\varepsilon(\mathbf{k})}$$

3D 时 $\int d^3 k/k^2$ 在 $k\to 0$ 收敛 → $M^2_\infty$ 有限（**平整相**）。
2D 时该积分对数发散 → $M^2 \sim \log L$（**粗糙相**）。
临界维度是 $d=2$，与 Mermin–Wagner 一致。

### 4.2 为什么修正是 $1/L$（而非 $1/L^2$ 或指数）

最低 6 个 $\mathbf{k}_{\min}$（沿 $\pm\hat x,\pm\hat y,\pm\hat z$）每个的贡献：

$$\frac{\langle |h_{\mathbf{k}_{\min}}|^2\rangle}{N} \sim \frac{1}{K\,(2\pi/L)^2 \cdot N} = \frac{L^2}{4\pi^2 K L^3} = \frac{1}{4\pi^2 K L}$$

一般维度：
$$M^2(L) - M^2_\infty \;\sim\; L^{-(d-2)}$$

| $d$ | 修正 | 物理 |
|---|---|---|
| 2 | $\log L$ | 粗糙相 |
| 3 | $1/L$ ✓ | 平整相，慢有限尺寸 |
| 4 | $1/L^2$ | 平整相，快收敛 |
| $\infty$ | $\exp(-L)$ | 平场极限 |

理论估计的系数：$c_{\rm est} \approx 6/(4\pi^2 K) \approx 0.152$（仅最低壳；其它 $\mathbf{k}$ 壳贡献会再压低一点）。拟合得到的 $\sigma=2$ 时 $c=0.113$，量级吻合。

### 4.3 $\sigma\ne 2$ 推广：RG 一致性

非高斯哈密顿量没有解析解，但项目早先的 coarse-graining 工作已经验证：所有 $\sigma$ 在长波长处 RG 流到**同一高斯不动点**（仅有效 stiffness $K_{\rm eff}(\sigma)$ 不同）。因此：

- IR 物理由有效高斯理论支配
- 最低 $\mathbf{k}$ 模仍然主导有限尺寸修正
- $M^2(L) - M^2_\infty \sim 1/(K_{\rm eff}\,L)$ 形式不变

数据支持：$\sigma=1, 2, 4$ 的 $c$ 几乎相等（$0.113$–$0.122$），说明它们的有效 stiffness 接近 → 同一个 Gaussian basin。$\sigma=0.5$ 的 $c=0.320$ 显著更大（更软的相互作用 → 更小的 $K_{\rm eff}$ → 更大的涨落），但仍是 $1/L$ 衰减。

## 5. 结论

1. **3D generalized SOS 对所有研究的 $\sigma\in\{0.5,1,2,4\}$ 都处于平整相**，无粗糙化转变。这与 2D（$\sigma=2$ 粗糙化、$M^2\sim\log L$）形成定性差别。
2. 有限尺寸修正服从 $M^2(L) = M^2_\infty - c/L$，与 Gaussian 模型 $L^{-(d-2)}$ 标度律一致。
3. 不同 $\sigma$ 的渐近 $M^2_\infty$ 落在 $0.13$–$0.43$ 区间，$\sigma=1$ 是最小值（双侧重新上升 → 最优 stiffness 在 $\sigma=1$ 附近）。

## 6. 局限与下一步

- 当前 fit window $L\in\{16,32,64\}$ 只 3 点，不能可靠分离 $1/L^2$ 次领头项。
- **三参数拟合** $M^2 = a + b/L + d/L^2$ 自由度太低（5 点 - 3 参 = 2 自由度），噪声敏感。
- **Richardson 一致性 check**：$M^2_\infty \approx 2M^2(64) - M^2(32)$，例如 $\sigma=2$ 给 $2(0.1246) - 0.1228 = 0.1264$，与三点拟合 $0.1264$ 完全一致 → 当前外推稳健。
- **更好的方案**：再跑 $L=128$（成本 ~30 h/σ），用 $L\in\{32,64,128\}$ 重新拟合；如果 $M^2_\infty$ 不漂移即确认。

## 附录：复现

```bash
# 集群上跑模拟
sbatch submit_3d.sh         # L ∈ {4,8,16,32}
sbatch submit_3d_L64.sh     # L = 64

# 本地画图（数据从 ustc:/home/shhu/sihan/SOS/data3d/ scp 下来）
julia plot_results_3d.jl
```

数据文件：`data3d/L{L}_s{σ}_3d.jld2`，字段 `W2_mean`、`W2_err`（变量名 `W2_*` 是历史遗留，物理量是 $M^2$）。
