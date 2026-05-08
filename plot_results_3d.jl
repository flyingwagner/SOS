using JLD2
using PyCall
using PyPlot
using Statistics

plt = pyimport("matplotlib.pyplot")
try
    pyimport("scienceplots")
    plt.style.use(["science"])
catch
    println("scienceplots 未安装，使用默认样式")
end
ticker = pyimport("matplotlib.ticker")

const DATA_DIR = "data3d"
const FIG_DIR  = "figures"
mkpath(FIG_DIR)

L_list     = [4, 8, 16, 32, 64]
sigma_list = [0.5, 1.0, 2.0, 4.0]

# ----- 读取 -----
function load_M2(L::Int, sigma::Float64)
    f = joinpath(DATA_DIR, "L$(L)_s$(sigma)_3d.jld2")
    d = load(f)
    return d["W2_mean"], d["W2_err"]   # 字段名 W2_* 是历史遗留，物理量是 M²
end

# 仅用 L ≥ 16 的点做 1/L 线性外推（小 L 偏离渐近过大）
fit_L_min = 16

# ----- 双面板图：M²(L) 与 M²(1/L) 外推 -----
fig, (ax1, ax2) = subplots(1, 2, figsize=(11, 4.5))

# ----- (a) M² vs L, log-x -----
ax1.set_xlabel("System size \$L\$")
ax1.set_ylabel("\$M^2\$")
ax1.set_xscale("log")

# ----- (b) M² vs 1/L 线性外推 -----
ax2.set_xlabel("\$1/L\$")
ax2.set_ylabel("\$M^2\$")

println("\n--- 3D 平整相 1/L 外推 (M² = M²_∞ - c/L)，仅使用 L ≥ $fit_L_min ---")

colors = ["C0","C1","C2","C3"]
for (i, sigma) in enumerate(sigma_list)
    M2_vals = [load_M2(L, sigma)[1] for L in L_list]
    M2_errs = [load_M2(L, sigma)[2] for L in L_list]

    # 左图
    ax1.errorbar(L_list, M2_vals, yerr=M2_errs, marker="o", color=colors[i],
                 label="\$\\sigma=$(sigma)\$", capsize=3, lw=1.2)

    # 右图：x = 1/L
    invL_all = 1.0 ./ L_list
    ax2.errorbar(invL_all, M2_vals, yerr=M2_errs, marker="o", color=colors[i],
                 label="\$\\sigma=$(sigma)\$", capsize=3, lw=0)

    # 线性拟合：仅 L ≥ fit_L_min
    mask = L_list .>= fit_L_min
    x_fit = 1.0 ./ L_list[mask]
    y_fit = M2_vals[mask]

    # 直接用 polyfit / linregress
    n = length(x_fit)
    xbar = mean(x_fit); ybar = mean(y_fit)
    slope = sum((x_fit .- xbar) .* (y_fit .- ybar)) / sum((x_fit .- xbar).^2)
    intercept = ybar - slope * xbar
    M2_inf = intercept

    # 拟合线延伸到 1/L = 0
    x_line = collect(range(0.0, maximum(invL_all)*1.05, length=50))
    y_line = intercept .+ slope .* x_line
    ax2.plot(x_line, y_line, "--", color=colors[i], lw=0.9, alpha=0.7)
    ax2.plot([0.0], [M2_inf], marker="*", color=colors[i], markersize=11,
             markeredgecolor="black", markeredgewidth=0.5, zorder=10)

    println("  σ=$sigma : M²_∞ = $(round(M2_inf, digits=4))   (slope c = $(round(-slope, digits=4)))")
end

# 左图美化
ax1.set_xticks(L_list)
ax1.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
ax1.set_xlim(3.5, maximum(L_list)*1.2)
ax1.legend(loc="best")
ax1.set_title("(a) \$M^2(L)\$")

# 右图：x 从 0 开始
ax2.set_xlim(-0.01, maximum(1.0 ./ L_list) * 1.05)
ax2.legend(loc="best")
ax2.set_title("(b) \$M^2\$ vs \$1/L\$, extrapolation to \$L\\to\\infty\$")

fig.tight_layout()
out = joinpath(FIG_DIR, "M2_scaling_3d.pdf")
fig.savefig(out)
println("\n已保存 $out")
