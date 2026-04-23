using JLD2
using PyCall
using PyPlot
using Random
using Statistics

plt = pyimport("matplotlib.pyplot")
pyimport("scienceplots")
plt.style.use(["science"])
optimize = pyimport("scipy.optimize")
scipy_stats = pyimport("scipy.stats")

const DATA_DIR   = "data_extra"
const FIG_DIR    = "figures"
const HIST_DIR   = joinpath(FIG_DIR, "histogram")
const SIGMA_LIST = [0.5, 1.0, 2.0, 4.0]
const L_LIST     = [64, 256, 512]

mkpath(HIST_DIR)

# ==========================================
# 辅助：提取 k 空间散点 + 大 k 区域 bin 平均
# ==========================================
"""
将 (k_abs, hk2) 分成两部分:
- k < k_cut: 全部保留
- k >= k_cut: 随机抽样 n_sample 个点
返回 (k_lo, hk2_lo, k_hi, hk2_hi)
"""
function split_and_subsample(k_abs, hk2_vals; k_cut=0.5, n_sample=200)
    lo = k_abs .< k_cut
    k_lo   = k_abs[lo]
    hk2_lo = hk2_vals[lo]

    hi = findall(k_abs .>= k_cut)
    if length(hi) > n_sample
        idx = hi[randperm(length(hi))[1:n_sample]]
    else
        idx = hi
    end
    sort!(idx)  # 保持 k 顺序
    return k_lo, hk2_lo, k_abs[idx], hk2_vals[idx]
end

"""从 hk_power_avg 矩阵提取 (k_abs, hk2_vals) 散点（全部 2D k 点）"""
function extract_k_hk2(hk2::Matrix{Float64}, L::Int)
    kx = [(i <= L÷2 ? i-1 : i-1-L) * (2π/L) for i in 1:L]
    ky = [(j <= L÷2 ? j-1 : j-1-L) * (2π/L) for j in 1:L]
    k_abs = Float64[]
    hk2_vals = Float64[]
    for j in 1:L, i in 1:L
        kmag = sqrt(kx[i]^2 + ky[j]^2)
        kmag == 0 && continue
        push!(k_abs, kmag)
        push!(hk2_vals, hk2[i, j])
    end
    return k_abs, hk2_vals
end

"""从 hk_power_avg 矩阵只提取 x 方向 (ky=0) 的 (kx, hk2) 数据"""
function extract_k_hk2_xdir(hk2::Matrix{Float64}, L::Int)
    kx_vals = Float64[]
    hk2_vals = Float64[]
    for i in 2:L  # 跳过零模 (i=1, j=1)
        kx = (i <= L÷2 ? i-1 : i-1-L) * (2π/L)
        push!(kx_vals, abs(kx))
        push!(hk2_vals, hk2[i, 1])  # j=1 对应 ky=0
    end
    return kx_vals, hk2_vals
end

"""
对 k < k_fit_max 的数据做 log-log OLS 线性回归: log(hk2) = log(A) - γ·log(k)
返回 NamedTuple (A, A_err, γ, γ_err, R2)
"""
function fit_power_law(k_abs, hk2_vals; k_fit_max=0.5)
    mask = k_abs .< k_fit_max
    lk = log.(k_abs[mask])
    lh = log.(hk2_vals[mask])
    result = scipy_stats.linregress(lk, lh)
    slope     = result[1]
    intercept = result[2]
    r         = result[3]
    slope_err = result[5]
    n = length(lk)
    s2 = sum((lh .- (slope .* lk .+ intercept)).^2) / (n - 2)
    inter_err = sqrt(s2 * (1.0/n + mean(lk)^2 / sum((lk .- mean(lk)).^2)))
    A     = exp(intercept)
    A_err = A * inter_err
    γ     = -slope
    γ_err = slope_err
    return (A=A, A_err=A_err, γ=γ, γ_err=γ_err, R2=r^2)
end

# ==========================================
# 拟合模型: P(x) = A * exp(-B * x^α)
# ==========================================
function fit_exp_model(x, y)
    # log P = log A - B * x^α
    # 参数: [A, B, α]
    model(x, A, B, α) = A .* exp.(-B .* x .^ α)
    p0 = [maximum(y), 1.0, 1.0]
    popt, _ = optimize.curve_fit(model, x, y, p0=p0,
                                  bounds=([0.0, 0.0, 0.1], [Inf, Inf, 10.0]),
                                  maxfev=10000)
    return popt
end

# ==========================================
# 用 L512 数据拟合 histogram，每个 sigma 一张图
# 自由 α，只拟合核心区域（P > P_max * 10^{-DROP}），确保 |Δh|→0 处拟合准确
# 画线性化图：log P vs |Δh|^α → 核心区域为直线
# ==========================================
const FIT_L = 512
const FIT_DROP_DECADES = 2   # 只拟合 P > P_max * 10^{-DROP_DECADES} 的核心区域
println("\n--- Histogram 拟合结果 (L=$FIT_L): P(|Δh|) = A·exp(-B·|Δh|^α)  (自由α, 核心拟合) ---")

for sigma in SIGMA_LIST
    fname = joinpath(DATA_DIR, "L$(FIT_L)_s$(sigma)_extra.jld2")
    if !isfile(fname)
        println("  跳过: $fname 不存在")
        continue
    end

    d = load(fname)
    bin_centers = d["bin_centers"]
    hist_counts = d["hist_counts"]
    bw = d["BIN_WIDTH"]
    prob = hist_counts ./ (sum(hist_counts) * bw)

    # 核心区域拟合：只用 P > P_max * 10^{-DROP_DECADES}
    P_max = maximum(prob)
    P_cut = P_max * 10.0^(-FIT_DROP_DECADES)
    mask_fit = prob .> P_cut
    x_fit = bin_centers[mask_fit]
    y_fit = prob[mask_fit]

    # 全部有效数据用于画图
    mask_plot = prob .> 1e-10

    try
        popt = fit_exp_model(x_fit, y_fit)
        A, B, α = popt
        println("  σ=$sigma : A=$(round(A,digits=4)), B=$(round(B,digits=4)), α=$(round(α,digits=4))")
        println("    核心拟合范围: |Δh| ∈ [$(round(minimum(x_fit),digits=3)), $(round(maximum(x_fit),digits=3))], P_cut=$(round(P_cut,digits=4))")

        # 线性化图: x = |Δh|^α, y = P (log scale)
        fig, ax = subplots(figsize=(6, 5))
        ax.set_title("\$L=$FIT_L, \\; \\sigma=$sigma, \\; \\alpha=$(round(α,digits=2))\$")
        ax.set_xlabel("\$|\\Delta h|^{$(round(α,digits=2))}\$")
        ax.set_ylabel("\$P(|\\Delta h|)\$")
        ax.set_yscale("log")

        # 数据点
        x_plot = bin_centers[mask_plot] .^ α
        y_plot = prob[mask_plot]
        ax.plot(x_plot, y_plot, "o", markersize=2, alpha=0.6, label="MC data")

        # 拟合直线 (在 |Δh|^α 空间中为直线)
        x_line = range(0, stop=maximum(x_plot), length=200)
        y_line = A .* exp.(-B .* collect(x_line))
        ax.plot(x_line, y_line, "-", lw=1.5, color="red",
                label="\$$(round(A,digits=2))\\,e^{-$(round(B,digits=2))\\,|\\Delta h|^{$(round(α,digits=2))}}\$")

        # 标记拟合范围
        x_cut = maximum(x_fit) ^ α
        ax.axvline(x_cut, color="gray", ls=":", lw=0.8, alpha=0.5, label="fit cutoff")

        ax.legend()
        fig.tight_layout()
        outpath = joinpath(HIST_DIR, "hist_L$(FIT_L)_s$(sigma).pdf")
        fig.savefig(outpath)
        plt.close(fig)
        println("  已保存 $outpath")
    catch e
        println("  σ=$sigma : 拟合失败 — $e")
    end
end

# ==========================================
# 汇总图：每个 L 一张，所有 sigma 叠加
# ==========================================
for L in L_LIST
    fig1, ax1 = subplots(figsize=(6, 5))
    ax1.set_title("Height Difference Distribution (\$L=$L\$)")
    ax1.set_xlabel("\$|\\Delta h|\$")
    ax1.set_ylabel("\$P(|\\Delta h|)\$")
    ax1.set_yscale("log")

    for sigma in SIGMA_LIST
        fname = joinpath(DATA_DIR, "L$(L)_s$(sigma)_extra.jld2")
        isfile(fname) || continue
        d = load(fname)
        prob = d["hist_counts"] ./ (sum(d["hist_counts"]) * d["BIN_WIDTH"])
        ax1.plot(d["bin_centers"], prob, label="\$\\sigma=$(sigma)\$")
    end

    ax1.legend()
    fig1.tight_layout()
    outpath = joinpath(FIG_DIR, "height_diff_histogram_L$(L).pdf")
    fig1.savefig(outpath)
    plt.close(fig1)
    println("\n已保存 $outpath")
end

# ==========================================
# 结构因子 ⟨|h_k|²⟩ vs k_x：每个 L 一张（只用 x 方向 ky=0）
# ==========================================
for Lval in L_LIST
    fig2, ax2 = subplots(figsize=(6, 5))
    ax2.set_title("Structure Factor (\$L=$Lval\$)")
    ax2.set_xlabel("\$|k_x|\$")
    ax2.set_ylabel("\$\\langle |h_{k_x}|^2 \\rangle\$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    println("\n--- 结构因子幂律拟合 (L=$Lval, x方向 ky=0, k < 0.5): ⟨|h_k|²⟩ ~ A/|k|^γ ---")
    for (ci, sigma) in enumerate(SIGMA_LIST)
        fname = joinpath(DATA_DIR, "L$(Lval)_s$(sigma)_extra.jld2")
        isfile(fname) || continue
        d = load(fname)
        L = d["L"]
        hk2 = d["hk_power_avg"]

        # 只提取 x 方向 (ky=0)
        k_abs, hk2_vals = extract_k_hk2_xdir(hk2, L)

        color = "C$(ci-1)"
        # 如果有误差数据则画 errorbar
        if haskey(d, "hk_power_err")
            hk_err = d["hk_power_err"]
            _, err_vals = extract_k_hk2_xdir(hk_err, L)
            ax2.errorbar(k_abs, hk2_vals, yerr=err_vals, fmt="o", markersize=4,
                         alpha=0.7, color=color, capsize=2,
                         label="\$\\sigma=$(sigma)\$")
        else
            ax2.plot(k_abs, hk2_vals, "o", markersize=4, alpha=0.7, color=color,
                     label="\$\\sigma=$(sigma)\$")
        end

        # 幂律拟合（只用小 k 数据）
        f = fit_power_law(k_abs, hk2_vals)
        println("  σ=$sigma : γ=$(round(f.γ,digits=4))±$(round(f.γ_err,digits=4)), A=$(round(f.A,digits=4))±$(round(f.A_err,digits=4)), R²=$(round(f.R2,digits=4))")
        k_fit = range(minimum(k_abs), stop=maximum(k_abs), length=200)
        ax2.plot(k_fit, f.A ./ k_fit .^ f.γ, "--", color=color, lw=1.2, alpha=0.8,
                 label="\$$(round(f.A,digits=2))/k^{$(round(f.γ,digits=2))}\$")
    end

    ax2.legend(fontsize=7)
    fig2.tight_layout()
    outpath = joinpath(FIG_DIR, "structure_factor_full_L$(Lval).pdf")
    fig2.savefig(outpath)
    plt.close(fig2)
    println("已保存 $outpath")
end

# ==========================================
# 固定 sigma，不同 L 对比 histogram
# ==========================================
for sigma in SIGMA_LIST
    fig, ax = subplots(figsize=(6, 5))
    ax.set_title("Height Difference Distribution (\$\\sigma=$sigma\$)")
    ax.set_xlabel("\$|\\Delta h|\$")
    ax.set_ylabel("\$P(|\\Delta h|)\$")
    ax.set_yscale("log")

    for L in L_LIST
        fname = joinpath(DATA_DIR, "L$(L)_s$(sigma)_extra.jld2")
        isfile(fname) || continue
        d = load(fname)
        prob = d["hist_counts"] ./ (sum(d["hist_counts"]) * d["BIN_WIDTH"])
        ax.plot(d["bin_centers"], prob, label="\$L=$L\$")
    end

    ax.legend()
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "height_diff_histogram_s$(sigma).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("已保存 $outpath")
end

# ==========================================
# 固定 sigma，不同 L 对比结构因子（径向 bin 平均）
# ==========================================
for sigma in SIGMA_LIST
    fig, ax = subplots(figsize=(6, 5))
    ax.set_title("Structure Factor (\$\\sigma=$sigma\$)")
    ax.set_xlabel("\$|k|\$")
    ax.set_ylabel("\$\\langle |h_k|^2 \\rangle\$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    # 反序画：大 L 先画（底层），小 L 后画（顶层）
    L_rev = reverse(L_LIST)
    nL = length(L_LIST)
    for (ri, Lval) in enumerate(L_rev)
        ci = nL - ri  # 保持颜色顺序: L_LIST[1]→C0, L_LIST[2]→C1, ...
        fname = joinpath(DATA_DIR, "L$(Lval)_s$(sigma)_extra.jld2")
        isfile(fname) || continue
        d = load(fname)
        L = d["L"]
        hk2 = d["hk_power_avg"]

        k_abs, hk2_vals = extract_k_hk2(hk2, L)

        color = "C$ci"
        ms = ri == nL ? 5 : 3  # 最小 L 的点稍大
        k_lo, hk2_lo, k_hi, hk2_hi = split_and_subsample(k_abs, hk2_vals)
        ax.scatter(k_lo, hk2_lo, s=ms*2, alpha=0.7, color=color, label="\$L=$Lval\$", zorder=ri+2)
        ax.plot(k_hi, hk2_hi, "o", markersize=ms-1, alpha=0.2, color=color, zorder=ri+2)
    end

    # 用最大 L 的数据做拟合
    fname_fit = joinpath(DATA_DIR, "L$(L_LIST[end])_s$(sigma)_extra.jld2")
    if isfile(fname_fit)
        d = load(fname_fit)
        k_abs, hk2_vals = extract_k_hk2(d["hk_power_avg"], d["L"])
        f = fit_power_law(k_abs, hk2_vals)
        println("  σ=$sigma (fit L=$(L_LIST[end])): γ=$(round(f.γ,digits=4))±$(round(f.γ_err,digits=4)), A=$(round(f.A,digits=4))±$(round(f.A_err,digits=4)), R²=$(round(f.R2,digits=4))")
        k_fit = range(minimum(k_abs), stop=π, length=200)
        ax.plot(k_fit, f.A ./ k_fit .^ f.γ, "k--", lw=1.2, alpha=0.8,
                label="\$$(round(f.A,digits=2))/k^{$(round(f.γ,digits=2))}\$", zorder=nL+3)
    end

    ax.legend(fontsize=7)
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "structure_factor_s$(sigma).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("已保存 $outpath")
end
