using JLD2
using PyCall
using PyPlot

plt = pyimport("matplotlib.pyplot")
pyimport("scienceplots")
plt.style.use(["science"])
ticker = pyimport("matplotlib.ticker")

const DATA_FILE = "sos_data.jld2"
const FIG_DIR   = "figures"
mkpath(FIG_DIR)

# ==========================================
# 读取数据
# ==========================================
data = load(DATA_FILE, "results")
meta = data["_meta"]

L_list     = meta["L_list"]
sigma_list = meta["sigma_list"]

# ==========================================
# 辅助：从 data 中取均值和误差
# ==========================================
get_W2(L, sigma)   = data["L$(L)_s$(sigma)"]["W2_mean"],   data["L$(L)_s$(sigma)"]["W2_err"]
get_chi(L, sigma)  = data["L$(L)_s$(sigma)"]["chi_mean"],  data["L$(L)_s$(sigma)"]["chi_err"]
get_chip(L, sigma) = data["L$(L)_s$(sigma)"]["chip_mean"], data["L$(L)_s$(sigma)"]["chip_err"]

scipy_stats = pyimport("scipy.stats")

# ==========================================
# 图1：W²(L) 有限尺寸标度
# ==========================================
fig1, ax1 = subplots(figsize=(6, 5))
ax1.set_title("Surface Roughness \$W^2(L)\$")
ax1.set_xlabel("System Size \$L\$")
ax1.set_ylabel("\$W^2\$")
ax1.set_xscale("log")

println("\n--- FSS 拟合结果 (W² ~ a·ln(L) + b) ---")
for sigma in sigma_list
    W2_vals = [get_W2(L, sigma)[1] for L in L_list]
    W2_errs = [get_W2(L, sigma)[2] for L in L_list]
    ax1.errorbar(L_list, W2_vals, yerr=W2_errs, marker="o", label="\$\\sigma=$(sigma)\$", capsize=4)

    slope, intercept, r, _, stderr = scipy_stats.linregress(log.(L_list), W2_vals)
    println("  σ=$sigma : W² ~ $(round(slope,digits=4))·ln(L) + $(round(intercept,digits=4))  (R²=$(round(r^2,digits=4)))")
end
ax1.set_xticks(sort(L_list))
ax1.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
ax1.set_xlim(3, maximum(L_list) * 1.2)
ax1.legend()
fig1.tight_layout()
fig1.savefig(joinpath(FIG_DIR, "W2_scaling.pdf"))
println("已保存 $(FIG_DIR)/W2_scaling.pdf")

# ==========================================
# 图2：χ_kmin(L) 随尺寸变化
# ==========================================
fig2, ax2 = subplots(figsize=(6, 5))
ax2.set_title("Structure Factor \$\\chi_{k_{\\min}}(L)\$")
ax2.set_xlabel("System Size \$L\$")
ax2.set_ylabel("\$\\chi_{k_{\\min}}\$")
ax2.set_xscale("log")
ax2.set_yscale("log")

println("\n--- χ_kmin 幂律拟合 (χ ~ L^α) ---")
for sigma in sigma_list
    chi_vals = [get_chi(L, sigma)[1] for L in L_list]
    chi_errs = [get_chi(L, sigma)[2] for L in L_list]
    ax2.errorbar(L_list, chi_vals, yerr=chi_errs, marker="s", label="\$\\sigma=$(sigma)\$", capsize=4)

    slope, intercept, r, _, _ = scipy_stats.linregress(log.(L_list), log.(chi_vals))
    println("  σ=$sigma : χ ~ L^$(round(slope,digits=4))  (R²=$(round(r^2,digits=4)))")
end
ax2.set_xticks(sort(L_list))
ax2.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
ax2.set_xlim(3, maximum(L_list) * 1.2)
ax2.legend()
fig2.tight_layout()
fig2.savefig(joinpath(FIG_DIR, "chi_scaling.pdf"))
println("已保存 $(FIG_DIR)/chi_scaling.pdf")

# ==========================================
# 图3：χ'_kmin(L) 随尺寸变化
# ==========================================
fig3, ax3 = subplots(figsize=(6, 5))
ax3.set_title("Connected Structure Factor \$\\chi'_{k_{\\min}}(L)\$")
ax3.set_xlabel("System Size \$L\$")
ax3.set_ylabel("\$\\chi'_{k_{\\min}}\$")
ax3.set_xscale("log")
ax3.set_yscale("log")

println("\n--- χ'_kmin 幂律拟合 (χ' ~ L^α) ---")
for sigma in sigma_list
    chip_vals = [get_chip(L, sigma)[1] for L in L_list]
    chip_errs = [get_chip(L, sigma)[2] for L in L_list]
    ax3.errorbar(L_list, chip_vals, yerr=chip_errs, marker="^", label="\$\\sigma=$(sigma)\$", capsize=4)

    slope, intercept, r, _, _ = scipy_stats.linregress(log.(L_list), log.(chip_vals))
    println("  σ=$sigma : χ' ~ L^$(round(slope,digits=4))  (R²=$(round(r^2,digits=4)))")
end
ax3.set_xticks(sort(L_list))
ax3.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
ax3.set_xlim(3, maximum(L_list) * 1.2)
ax3.legend()
fig3.tight_layout()
fig3.savefig(joinpath(FIG_DIR, "chip_scaling.pdf"))
println("已保存 $(FIG_DIR)/chip_scaling.pdf")
