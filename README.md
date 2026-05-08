# Generalized Continuous SOS Model — Monte Carlo Study

二维连续 Solid-on-Solid (SOS) 模型的蒙特卡洛模拟、有限尺寸标度分析、结构因子、以及实空间块平均（coarse-graining）下向高斯不动点的 RG 流动。

Hamiltonian:

$$\mathcal{H} = K \sum_{\langle i,j \rangle} |h_i - h_j|^{\sigma}, \qquad \sum_i h_i = 0$$

$L \in \{4,8,16,32,64,128,256,512\}$，$\sigma \in \{0.5, 1.0, 2.0, 4.0\}$。Coarse-grain 分析在 $L=1024$ 上做 $2\times2$ 块平均的多级 RG 流。

---

## Codebase

### Core library

| File | 作用 |
|---|---|
| `GeneralizedSOS.jl` | 哈密顿量、近邻枚举、Metropolis 更新、零模约束 (`Σh_i = 0`)。 |
| `SOSCore.jl`        | 2D 方格子 MC 主循环、热化、测量间隔、自相关控制。 |
| `SOSCore3D.jl`      | 3D 立方格子版本（6 近邻、嵌入反射 SW 集团、PBC、零模约束）。 |
| `Measure.jl`        | 2D 单次构型观测量：$M^2$、$G(\mathbf{r})$、$h_{\mathbf{k}}$、能量密度等。 |
| `Measure3D.jl`      | 3D 单次构型观测量：$M^2$、$G(r)$（沿 x）、$\chi_{\mathbf{k}_{\min}}$、$|\Delta h|$ 直方图、3D FFT 全谱。 |
| `BlockStats.jl`     | 蒙特卡洛测量的分块误差估计（block bootstrap / jackknife）。 |
| `CoarseGrain.jl`    | $2\times2$ 实空间块平均，用于 RG 流分析。 |

### Drivers

| File | 作用 |
|---|---|
| `run.jl`                 | 2D 主模拟入口：扫描 $(L,\sigma)$，输出到 `data/` 下 `L{L}_s{σ}.jld2`。 |
| `run3d.jl`               | 3D 主模拟入口：$L^3$ 立方格子，输出到 `data3d/` 下 `L{L}_s{σ}_3d.jld2`。命令行同 `run.jl`。 |
| `run_measure_extra.jl`   | 在已有平衡态上追加测量：$h_{\mathbf{k}}$ 的 FFT 幂谱及其 sum-of-squares 误差，输出到 `data_extra/`。命令行接受 `L σ`。 |
| `run_coarse_grain.jl`    | $L=1024$ 的 coarse-graining 流水线，逐级块平均并累计各阶矩与直方图，输出到 `data_cg/`。 |
| `merge_data.jl`          | 合并多个随机种子/分片的数据文件。 |
| `bench.jl`, `benchmark_histogram.jl` | 性能基准。 |

### Analysis & plots

| File | 作用 |
|---|---|
| `fit_analysis.jl`         | $M^2(L)$、$\chi_{\mathbf{k}}(L)$ 等关于 $L$ 的标度拟合（对数-对数回归）。 |
| `plot_results.jl`         | 主结果图：$M^2$-$L$、$\chi$、$\chi'$ 的有限尺寸标度。 |
| `plot_histogram_chik.jl`  | 高度差 $|\Delta h|$ 分布直方图、伸展指数拟合 $P(|\Delta h|) = A e^{-B |\Delta h|^{\alpha}}$；结构因子 $S(k_x, k_y=0)$ 及 x-方向误差棒。 |
| `plot_coarse_grain.jl`    | Coarse-grain 级数 $b = 1, 2, 4, \ldots$ 下 $|\Delta h|^2/w^2$ 的数据坍缩到 half-normal。 |

### SLURM (ihicnormal cluster)

`submit_extra.sh`, `submit_cg.sh`, `run_extra_all.sh`, `test_job.sh`：批量提交脚本，均已按 `--partition=ihicnormal` 配置并使用绝对路径。

### Report

`report/report.tex` — LaTeX 源。`report/report.pdf` / `report/report_SOS.pdf` — 编译产物。`report/review.md` — review 意见记录。

### Documentation (`doc/`)

| 文件 | 内容 |
|---|---|
| `doc/simulation.md`         | 2D 模拟流程、观测量定义、运行参数说明。 |
| `doc/further_simulations.md`| 待补跑的 $(L,\sigma)$ 组合列表（`hk_power_err` 字段补齐）。 |
| `doc/analysis_3d.md`        | 3D 平整相 $M^2(L)$ 标度分析、$1/L$ 外推、Gaussian 不动点 RG 一致性论证。 |

---

## 数据目录（不纳入版本控制）

实际运行产出的 `.jld2` 数据 (~100MB+) 以及 `figures/` 下 PDF 都被 `.gitignore` 排除。需要时按下表重跑：

| 目录 | 内容 | 生成命令 |
|---|---|---|
| `data/`       | 2D 主模拟：$M^2$、$G$、$\chi_{\mathbf{k}}$ for $(L,\sigma)$ | `julia run.jl L σ` |
| `data3d/`     | 3D 立方格子模拟（同样观测量）                         | `julia run3d.jl L σ` |
| `data_extra/` | FFT 幂谱及误差                                       | `julia run_measure_extra.jl L σ` |
| `data_cg/`    | $L=1024$ 的 CG 级数数据                               | `julia run_coarse_grain.jl` |
| `figures/`    | 所有 plot 脚本输出的 PDF                              | 对应 `plot_*.jl` |

`doc/further_simulations.md` 列出尚需重跑以补齐 `hk_power_err` 字段的 $(L,\sigma)$ 组合。

---

## 项目执行情况

**已完成：**

1. 主模拟循环和核心观测量（$M^2$、$G(\mathbf{r})$、$\chi_{\mathbf{k}}$、能量）实现并跑通 $(L,\sigma) \in \{4,\dots,512\} \times \{0.5,1.0,2.0,4.0\}$。
2. 有限尺寸标度拟合：$M^2(L) \sim \log L$ 的粗糙相行为在 $\sigma = 2$ 处验证，其它 $\sigma$ 下提取有效指数。
3. 高度差分布：用自由 $\alpha$ 的伸展指数 $P(|\Delta h|) = A e^{-B|\Delta h|^\alpha}$ 拟合，核心区 $P > P_{\max}\cdot 10^{-2}$ 截断，$|\Delta h| \to 0$ 贴合良好。拟合参数：

   | σ   | A      | B      | α      |
   |-----|--------|--------|--------|
   | 0.5 | 0.9150 | 1.0045 | 0.8059 |
   | 1.0 | 1.5636 | 1.6102 | 1.3252 |
   | 2.0 | 1.5958 | 2.0000 | 2.0000 |
   | 4.0 | 1.4678 | 2.0422 | 2.7172 |

   $\alpha \ne \sigma$（除 $\sigma=2$ 外）源于边际分布 $\ne$ 单键 Boltzmann 权重。
4. 结构因子改为沿 $k_x$ 方向（$k_y=0$）绘制，并在测量中累计 sum-of-squares 以给出误差棒（`run_measure_extra.jl`）。
5. Coarse-graining RG 流：$L=1024$ 各 $\sigma$ 下多级 $2\times2$ 块平均，重标度后 $|\Delta \bar h|^2/w^2$ 所有级数坍缩到 half-normal（高斯不动点）。
6. 报告 `report.tex` 已同步以上全部修改并重新编译。

**待办：**

- 重新运行 `run_measure_extra.jl` 覆盖全部 12 组 $(L,\sigma)$ 以使 `hk_power_err` 字段写入所有数据文件（见 `doc/further_simulations.md`）；目前只有部分文件有误差数据，因此 Fig 4 的误差棒在部分 $\sigma$ 上缺失。

---

## Quick start

```bash
# 主模拟（单组 L, σ 示例需手动编辑 run.jl 里的扫描列表）
julia --project=. run.jl

# 追加 FFT 测量
julia --project=. run_measure_extra.jl 512 2.0

# Coarse-grain
julia --project=. run_coarse_grain.jl

# 画图
julia --project=. plot_results.jl
julia --project=. plot_histogram_chik.jl
julia --project=. plot_coarse_grain.jl

# 编译报告
cd report && pdflatex report.tex
```

集群（ihicnormal 分区）：`sbatch submit_extra.sh` / `sbatch submit_cg.sh`。
