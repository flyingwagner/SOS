using JLD2
using PyCall
using PyPlot

plt = pyimport("matplotlib.pyplot")
pyimport("scienceplots")
plt.style.use(["science"])
ticker   = pyimport("matplotlib.ticker")
colors   = pyimport("matplotlib.colors")
cm       = pyimport("matplotlib.cm")
optimize = pyimport("scipy.optimize")

const DATA_DIR   = "data_cg"
const FIG_DIR    = joinpath("figures", "coarse_grain")
mkpath(FIG_DIR)

# ==========================================
# 参数：扫描 data_cg 目录下所有文件
# ==========================================
all_files = filter(f -> endswith(f, "_cg.jld2"), readdir(DATA_DIR))
if isempty(all_files)
    error("data_cg/ 目录下没有找到 *_cg.jld2 文件")
end

# 解析 (L, sigma) 列表
function parse_filename(fname)
    m = match(r"L(\d+)_s([\d.]+)_cg\.jld2", fname)
    m === nothing && return nothing
    return parse(Int, m[1]), parse(Float64, m[2])
end

parsed = filter(!isnothing, [parse_filename(f) for f in all_files])
L_set     = sort(unique([p[1] for p in parsed]))
sigma_set = sort(unique([p[2] for p in parsed]))

println("找到数据: L = $L_set, σ = $sigma_set")

# ==========================================
# 辅助：加载单个文件
# ==========================================
function load_cg(L, sigma)
    fname = joinpath(DATA_DIR, "L$(L)_s$(sigma)_cg.jld2")
    isfile(fname) || return nothing
    return load(fname)
end

# ==========================================
# 辅助：颜色映射（按层级 b）
# ==========================================
function level_colors(n)
    cmap = cm.viridis
    return [cmap(i / (n - 1)) for i in 0:(n-1)]
end

# ==========================================
# 图1：⟨(Δh)²⟩ vs block size b（RG flow）
# 每个 sigma 一张，不同 L 叠加
# ==========================================
println("\n--- 图1: ⟨(Δh)²⟩ vs b ---")
for sigma in sigma_set
    fig, ax = subplots(figsize=(5, 4))
    ax.set_title("\$\\langle (\\Delta h)^2 \\rangle\$ vs block size \$b\$ (\$\\sigma=$sigma\$)")
    ax.set_xlabel("Block size \$b\$")
    ax.set_ylabel("\$\\langle (\\Delta h)^2 \\rangle\$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    for L in L_set
        d = load_cg(L, sigma)
        d === nothing && continue
        bs = d["block_sizes"]
        dh2 = d["mean_dh2"]
        ax.plot(bs, dh2, "o-", markersize=4, label="\$L=$L\$")
    end

    ax.legend()
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "dh2_vs_b_s$(sigma).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("  已保存 $outpath")
end

# ==========================================
# 图2：excess kurtosis κ₄ vs block size b
# 每个 sigma 一张，不同 L 叠加
# ==========================================
println("\n--- 图2: κ₄ vs b ---")
for sigma in sigma_set
    fig, ax = subplots(figsize=(5, 4))
    ax.set_title("Excess kurtosis \$\\kappa_4\$ vs block size \$b\$ (\$\\sigma=$sigma\$)")
    ax.set_xlabel("Block size \$b\$")
    ax.set_ylabel("Excess kurtosis \$\\kappa_4\$")
    ax.set_xscale("log")
    ax.axhline(0.0, color="gray", lw=0.8, ls="--", label="Gaussian")

    for L in L_set
        d = load_cg(L, sigma)
        d === nothing && continue
        bs  = d["block_sizes"]
        kurt = d["kurtosis_excess"]
        ax.plot(bs, kurt, "o-", markersize=4, label="\$L=$L\$")
    end

    ax.legend()
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "kurtosis_vs_b_s$(sigma).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("  已保存 $outpath")
end

# ==========================================
# 图3：shell 平均结构因子 S(k) vs |k|，各 CG 层级
# 每个 (L, sigma) 一张
# ==========================================
println("\n--- 图3: S(k) vs |k| 各 CG 层级 ---")
for sigma in sigma_set, L in L_set
    d = load_cg(L, sigma)
    d === nothing && continue

    shell_S  = d["shell_S_avg"]      # Vector{Vector{Float64}}, length N_LEVELS_K
    shell_km = d["shell_k_means"]    # Vector{Vector{Float64}}
    shell_ct = d["shell_counts"]     # Vector{Vector{Int}}
    bs       = d["block_sizes"]
    n_lev    = length(shell_S)
    clrs     = level_colors(n_lev)

    fig, ax = subplots(figsize=(5, 4))
    ax.set_title("Structure factor \$S(k)\$ (\$L=$L,\\;\\sigma=$sigma\$)")
    ax.set_xlabel("\$|k|\$")
    ax.set_ylabel("\$S(k) = \\langle |\\tilde{h}_k|^2 \\rangle\$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    for lev in 1:n_lev
        km  = shell_km[lev]
        S   = shell_S[lev]
        cnt = shell_ct[lev]
        mask = cnt .> 0
        b = bs[lev]
        ax.plot(km[mask], S[mask], "o-", markersize=3, color=clrs[lev],
                label="\$b=$b\$")
    end

    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "Sk_levels_L$(L)_s$(sigma).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("  已保存 $outpath")
end

# ==========================================
# 图4：实空间高度差直方图随 CG 层级演化（x 轴归一化为 |Δh|²/w²）
# 每个 (L, sigma) 一张
# 归一化后高斯分布变为标准形式: ln P ∝ -|Δh|²/(2w²), 斜率 -1/2
# ==========================================
println("\n--- 图4: P(|Δh|) 各 CG 层级 (x = |Δh|²/w²) ---")
for sigma in sigma_set, L in L_set
    d = load_cg(L, sigma)
    d === nothing && continue

    hists      = d["hist"]          # Vector{Vector{Float64}}, length N_LEVELS
    bin_ctrs   = d["bin_centers"]   # Vector{Float64}
    bw         = d["BIN_WIDTH"]
    bs         = d["block_sizes"]
    mean_dh2   = d["mean_dh2"]
    n_lev      = length(hists)
    clrs       = level_colors(n_lev)

    fig, ax = subplots(figsize=(5, 4))
    ax.set_title("Height diff. distribution (\$L=$L,\\;\\sigma=$sigma\$)")
    ax.set_xlabel("\$|\\Delta \\bar{h}|^2 / w^2\$")
    ax.set_ylabel("\$P(|\\Delta \\bar{h}|)\$")
    ax.set_yscale("log")

    for lev in 1:n_lev
        h_raw = hists[lev]
        norm  = sum(h_raw) * bw
        norm == 0 && continue
        prob  = h_raw[1:end-1] ./ norm
        bc    = bin_ctrs[1:end-1]
        mask  = prob .> 1e-12
        b = bs[lev]
        w2 = mean_dh2[lev]   # w² = ⟨(Δh)²⟩
        # 归一化 x 轴: |Δh|² / w²
        x_rescaled = (bc[mask] .^ 2) ./ w2
        ax.plot(x_rescaled, prob[mask], "-", lw=1.0, color=clrs[lev],
                label="\$b=$b\$")
    end

    # 标准半正态参考线: P = sqrt(2/π)/w · exp(-|Δh|²/(2w²))
    # 在归一化坐标下 t = |Δh|²/w²: ln P = const - t/2
    # 用 b=64 的 w 作为参考振幅
    w64 = sqrt(mean_dh2[findfirst(==(64), bs)])
    t_ref = range(0.0, stop=12.0, length=200)
    y_ref = sqrt(2 / pi) / w64 .* exp.(-t_ref ./ 2)
    ax.plot(t_ref, y_ref, "k--", lw=1.4, alpha=0.6, label="half-normal")

    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "hist_levels_L$(L)_s$(sigma).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("  已保存 $outpath")
end

"""
Log-space 线性回归拟合: log P = a₀ + Σ aᵢ x^eᵢ
基函数 {x^σ, x, x²} 去重，只用 P > 1e-6 的点拟合。
返回 (exponents, coeffs) 其中 coeffs[1]=a₀, coeffs[2:end]=aᵢ
"""
function fit_logspace(x, prob, sigma)
    mask_fit = prob .> 1e-6
    xf = x[mask_fit]
    lp = log.(prob[mask_fit])
    exponents = sort(unique([sigma, 1.0, 2.0]))
    X = hcat(ones(length(xf)), [xf .^ e for e in exponents]...)
    beta = X \ lp
    return exponents, beta
end

"""用 logspace 拟合结果计算 P(x)"""
function eval_logspace(x, exponents, beta)
    lp = beta[1] .+ sum(beta[i+1] .* x .^ exponents[i] for i in eachindex(exponents))
    return exp.(lp)
end

"""格式化拟合公式为 LaTeX"""
function format_fit_label(exponents, beta)
    s = "\$$(round(exp(beta[1]),digits=2))"
    for (i, e) in enumerate(exponents)
        c = beta[i+1]
        c == 0 && continue
        sign = c < 0 ? "-" : "+"
        ac = round(abs(c), digits=2)
        if e == 1.0
            s *= "\\,e^{$(sign)$(ac)\\,x"
        elseif e == round(e)
            s *= "\\,e^{$(sign)$(ac)\\,x^{$(Int(e))}"
        else
            s *= "\\,e^{$(sign)$(ac)\\,x^{$(round(e,digits=1))}"
        end
    end
    # 合并成单个 exp
    parts = String[]
    for (i, e) in enumerate(exponents)
        c = beta[i+1]
        abs(c) < 1e-4 && continue
        sign = c < 0 ? "-" : (isempty(parts) ? "" : "+")
        ac = round(abs(c), digits=2)
        if e == 1.0
            push!(parts, "$(sign)$(ac)x")
        elseif e == round(e)
            push!(parts, "$(sign)$(ac)x^{$(Int(e))}")
        else
            push!(parts, "$(sign)$(ac)x^{$(round(e,digits=1))}")
        end
    end
    return "\$$(round(exp(beta[1]),digits=2))\\,e^{$(join(parts))}\$"
end

# ==========================================
# 图4b：每个 (sigma, b) 单独一张图，保存到 hist_real/
# ==========================================
const HIST_REAL_DIR = joinpath(FIG_DIR, "hist_real")
mkpath(HIST_REAL_DIR)

println("\n--- 图4b: 单独 histogram per (σ, b) ---")
for sigma in sigma_set, L in L_set
    d = load_cg(L, sigma)
    d === nothing && continue

    hists    = d["hist"]
    bin_ctrs = d["bin_centers"]
    bw       = d["BIN_WIDTH"]
    bs       = d["block_sizes"]
    n_lev    = length(hists)

    for lev in 1:n_lev
        h_raw = hists[lev]
        norm  = sum(h_raw) * bw
        norm == 0 && continue
        prob  = h_raw[1:end-1] ./ norm
        bc    = bin_ctrs[1:end-1]
        mask  = prob .> 1e-12
        b = bs[lev]

        fig, ax = subplots(figsize=(5, 4))
        ax.set_title("\$L=$L,\\;\\sigma=$sigma,\\;b=$b\$")
        ax.set_xlabel("\$|\\Delta h|\$")
        ax.set_ylabel("\$P(|\\Delta h|)\$")
        ax.set_yscale("log")

        ax.plot(bc[mask], prob[mask], "o", markersize=2, label="MC data")

        try
            exponents, beta = fit_logspace(bc[mask], prob[mask], sigma)
            x_dense = range(minimum(bc[mask]), stop=maximum(bc[mask]), length=300)
            y_dense = eval_logspace(collect(x_dense), exponents, beta)
            lbl = format_fit_label(exponents, beta)
            ax.plot(x_dense, y_dense, "-", lw=1.2, label=lbl)

            coef_str = join(["$(round(e,digits=1)):$(round(beta[i+1],digits=4))" for (i,e) in enumerate(exponents)], ", ")
            println("  σ=$sigma, b=$b: a0=$(round(beta[1],digits=4)), $coef_str")
        catch e
            println("  σ=$sigma, b=$b: 拟合失败 — $e")
        end

        ax.legend(fontsize=7)
        fig.tight_layout()
        outpath = joinpath(HIST_REAL_DIR, "hist_L$(L)_s$(sigma)_b$(b).pdf")
        fig.savefig(outpath)
        plt.close(fig)
    end
end
println("  hist_real/ 图像已保存")

# ==========================================
# 图4c：拟合系数 vs b（RG flow 到 Gaussian 不动点）
# 每个 sigma 一张
# ==========================================
println("\n--- 图4c: 拟合系数 vs b (RG flow) ---")
for sigma in sigma_set, L in L_set
    d = load_cg(L, sigma)
    d === nothing && continue

    hists    = d["hist"]
    bin_ctrs = d["bin_centers"]
    bw       = d["BIN_WIDTH"]
    bs       = d["block_sizes"]
    n_lev    = length(hists)

    # 收集所有层级的拟合系数
    exponents_all = Vector{Vector{Float64}}(undef, n_lev)
    coeffs_all    = Vector{Vector{Float64}}(undef, n_lev)

    for lev in 1:n_lev
        h_raw = hists[lev]
        norm  = sum(h_raw) * bw
        norm == 0 && continue
        prob  = h_raw[1:end-1] ./ norm
        bc    = bin_ctrs[1:end-1]
        mask  = prob .> 1e-12

        try
            exponents, beta = fit_logspace(bc[mask], prob[mask], sigma)
            exponents_all[lev] = exponents
            coeffs_all[lev]    = beta[2:end]  # 跳过 a0
        catch e
            exponents_all[lev] = Float64[]
            coeffs_all[lev]    = Float64[]
        end
    end

    # 画图：每个指数一条曲线
    fig, ax = subplots(figsize=(5, 4))
    ax.set_title("Fit coefficients vs block size (\$L=$L,\\;\\sigma=$sigma\$)")
    ax.set_xlabel("Block size \$b\$")
    ax.set_ylabel("Coefficient (absolute value)")
    ax.set_xscale("log")
    ax.set_yscale("log")

    # 找出所有出现过的指数
    all_exp = sort(unique(vcat(exponents_all...)))

    for exp in all_exp
        coef_vs_b = Float64[]
        bs_plot   = Int[]
        for lev in 1:n_lev
            if isempty(exponents_all[lev])
                continue
            end
            idx = findfirst(==(exp), exponents_all[lev])
            if idx !== nothing
                push!(coef_vs_b, abs(coeffs_all[lev][idx]))
                push!(bs_plot, bs[lev])
            end
        end
        if !isempty(coef_vs_b)
            if exp == 1.0
                lbl = "\$x\$"
            elseif exp == round(exp)
                lbl = "\$x^{$(Int(exp))}\$"
            else
                lbl = "\$x^{$(round(exp,digits=1))}\$"
            end
            ax.plot(bs_plot, coef_vs_b, "o-", markersize=4, label=lbl)
        end
    end

    ax.legend()
    ax.grid(true, alpha=0.3)
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "fit_coefs_vs_b_L$(L)_s$(sigma).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("  已保存 $outpath")
end

# ==========================================
# 图5：⟨(Δh)²⟩ vs b，固定 L，不同 sigma 对比
# 每个 L 一张
# ==========================================
println("\n--- 图5: ⟨(Δh)²⟩ vs b，固定 L ---")
for L in L_set
    fig, ax = subplots(figsize=(5, 4))
    ax.set_title("\$\\langle (\\Delta h)^2 \\rangle\$ vs \$b\$ (\$L=$L\$)")
    ax.set_xlabel("Block size \$b\$")
    ax.set_ylabel("\$\\langle (\\Delta h)^2 \\rangle\$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    for (ci, sigma) in enumerate(sigma_set)
        d = load_cg(L, sigma)
        d === nothing && continue
        bs  = d["block_sizes"]
        dh2 = d["mean_dh2"]
        ax.plot(bs, dh2, "o-", markersize=4, color="C$(ci-1)",
                label="\$\\sigma=$sigma\$")
    end

    ax.legend()
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "dh2_vs_b_L$(L).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("  已保存 $outpath")
end

# ==========================================
# 图6：κ₄ vs b，固定 L，不同 sigma 对比
# 每个 L 一张
# ==========================================
println("\n--- 图6: κ₄ vs b，固定 L ---")
for L in L_set
    fig, ax = subplots(figsize=(5, 4))
    ax.set_title("Excess kurtosis \$\\kappa_4\$ vs \$b\$ (\$L=$L\$)")
    ax.set_xlabel("Block size \$b\$")
    ax.set_ylabel("\$\\kappa_4\$")
    ax.set_xscale("log")
    ax.axhline(0.0, color="gray", lw=0.8, ls="--", label="Gaussian")

    for (ci, sigma) in enumerate(sigma_set)
        d = load_cg(L, sigma)
        d === nothing && continue
        bs   = d["block_sizes"]
        kurt = d["kurtosis_excess"]
        ax.plot(bs, kurt, "o-", markersize=4, color="C$(ci-1)",
                label="\$\\sigma=$sigma\$")
    end

    ax.legend()
    fig.tight_layout()
    outpath = joinpath(FIG_DIR, "kurtosis_vs_b_L$(L).pdf")
    fig.savefig(outpath)
    plt.close(fig)
    println("  已保存 $outpath")
end

println("\n全部图像已保存至 $FIG_DIR/")
