using JLD2
using PyCall
using PyPlot

plt = pyimport("matplotlib.pyplot")
pyimport("scienceplots")
plt.style.use(["science"])
optimize = pyimport("scipy.optimize")

const DATA_DIR   = "data_extra"
const FIG_DIR    = "figures"
const HIST_DIR   = joinpath(FIG_DIR, "histogram")
const SIGMA_LIST = [0.5, 1.0, 2.0, 4.0]
const L_LIST     = [64]

mkpath(HIST_DIR)

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
# 每个 (L, sigma) 单独画图 + 拟合
# ==========================================
println("\n--- Histogram 拟合结果: P(|Δh|) = A·exp(-B·|Δh|^α) ---")

for L in L_LIST
    for sigma in SIGMA_LIST
        fname = joinpath(DATA_DIR, "L$(L)_s$(sigma)_extra.jld2")
        if !isfile(fname)
            println("  跳过: $fname 不存在")
            continue
        end

        d = load(fname)
        bin_centers = d["bin_centers"]
        hist_counts = d["hist_counts"]
        bw = d["BIN_WIDTH"]
        prob = hist_counts ./ (sum(hist_counts) * bw)

        # 只拟合有效数据点 (prob > 0)
        mask = prob .> 1e-10
        x_fit = bin_centers[mask]
        y_fit = prob[mask]

        fig, ax = subplots(figsize=(6, 5))
        ax.set_title("\$L=$L, \\; \\sigma=$sigma\$")
        ax.set_xlabel("\$|\\Delta h|\$")
        ax.set_ylabel("\$P(|\\Delta h|)\$")

        ax.plot(bin_centers, prob, "o", markersize=3, label="MC data")

        try
            popt = fit_exp_model(x_fit, y_fit)
            A, B, α = popt
            println("  L=$L, σ=$sigma : A=$(round(A,digits=4)), B=$(round(B,digits=4)), α=$(round(α,digits=4))")

            x_dense = range(0, stop=maximum(bin_centers), length=500)
            y_dense = A .* exp.(-B .* x_dense .^ α)
            ax.plot(x_dense, y_dense, "-", lw=1.5,
                    label="\$$(round(A,digits=2)) e^{-$(round(B,digits=2))|\\Delta h|^{$(round(α,digits=2))}}\$")
        catch e
            println("  L=$L, σ=$sigma : 拟合失败 — $e")
        end

        ax.legend()
        fig.tight_layout()
        outpath = joinpath(HIST_DIR, "hist_L$(L)_s$(sigma).pdf")
        fig.savefig(outpath)
        plt.close(fig)
        println("  已保存 $outpath")
    end
end

# ==========================================
# 汇总图：所有 sigma 叠加
# ==========================================
fig1, ax1 = subplots(figsize=(6, 5))
ax1.set_title("Height Difference Distribution (\$L=$(L_LIST[1])\$)")
ax1.set_xlabel("\$|\\Delta h|\$")
ax1.set_ylabel("\$P(|\\Delta h|)\$")

for sigma in SIGMA_LIST
    fname = joinpath(DATA_DIR, "L$(L_LIST[1])_s$(sigma)_extra.jld2")
    isfile(fname) || continue
    d = load(fname)
    prob = d["hist_counts"] ./ (sum(d["hist_counts"]) * d["BIN_WIDTH"])
    ax1.plot(d["bin_centers"], prob, label="\$\\sigma=$(sigma)\$")
end

ax1.legend()
fig1.tight_layout()
fig1.savefig(joinpath(FIG_DIR, "height_diff_histogram.pdf"))
println("\n已保存 $(FIG_DIR)/height_diff_histogram.pdf")

# ==========================================
# 结构因子 ⟨|h_k|²⟩ vs |k|
# ==========================================
fig2, ax2 = subplots(figsize=(6, 5))
ax2.set_title("Structure Factor (\$L=$(L_LIST[1])\$)")
ax2.set_xlabel("\$|k|\$")
ax2.set_ylabel("\$\\langle |h_k|^2 \\rangle\$")
ax2.set_xscale("log")
ax2.set_yscale("log")

for sigma in SIGMA_LIST
    fname = joinpath(DATA_DIR, "L$(L_LIST[1])_s$(sigma)_extra.jld2")
    isfile(fname) || continue
    d = load(fname)
    L = d["L"]
    hk2 = d["hk_power_avg"]

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

    ax2.scatter(k_abs, hk2_vals, s=4, alpha=0.4, label="\$\\sigma=$(sigma)\$")
end

k_ref = range(0.1, stop=π, length=100)
ax2.plot(k_ref, 0.5 ./ k_ref.^2, "k--", lw=1.5, label="\$1/k^2\$")

ax2.legend()
fig2.tight_layout()
fig2.savefig(joinpath(FIG_DIR, "structure_factor_full.pdf"))
println("已保存 $(FIG_DIR)/structure_factor_full.pdf")
