include("GeneralizedSOS.jl")

using .GeneralizedSOS
using PyPlot
using PyCall
using Statistics

plt = pyimport("matplotlib.pyplot")
pyimport("scienceplots")
plt.style.use(["science"])

ticker = pyimport("matplotlib.ticker")

# ==========================================
# 模拟参数
# ==========================================
const CONV_SWEEPS    = 50000   # 任务1：收敛性检验总步数
const THERM_SWEEPS   = 50000   # 任务2：热化步数
const MEASURE_SWEEPS = 200000  # 任务2：测量步数

# ==========================================
# 任务1：L=16，sigma=1,2,3，画<m^2>随蒙卡步数变化
# 证明系统已收敛
# ==========================================

function run_convergence_check(; L::Int=16, sigma::Float64=2.0, K::Float64=1.0,
                                  total_sweeps::Int=CONV_SWEEPS, step_size::Float64=1.0)
    h = zeros(Float64, L, L)
    roughness_trace = Float64[]

    running_sum = 0.0
    for i in 1:total_sweeps
        GeneralizedSOS.metropolis_sweep!(h, L, K, sigma, step_size)
        if i % 5 == 0
            GeneralizedSOS.swendsen_wang_sweep!(h, L, K, sigma)
        end
        running_sum += GeneralizedSOS.measure_roughness(h)
        push!(roughness_trace, running_sum / i)
    end

    return roughness_trace
end

println("=== 任务1：收敛性检验 L=16 ===")
# sigma_list = [1.0, 2.0, 3.0]
sigma_list = [2.0]
total_sweeps = CONV_SWEEPS

fig1, ax1 = subplots(1, 1, figsize=(8, 5))

for sigma in sigma_list
    println("  运行 sigma=$sigma ...")
    trace = run_convergence_check(L=16, sigma=sigma, total_sweeps=total_sweeps)
    ax1.plot(1:total_sweeps, trace, label="\$\\sigma = $sigma\$", linewidth=0.8, alpha=0.85)
end

ax1.set_xlabel("Monte Carlo Sweeps", fontsize=13)
ax1.set_ylabel("\$\\langle m^2 \\rangle\$", fontsize=13)
ax1.set_title("Convergence Check: \$\\langle m^2 \\rangle\$ vs MC Sweeps, \$L=16\$", fontsize=13)
ax1.legend(fontsize=12)
ax1.grid(true, alpha=0.3)
fig1.tight_layout()
fig1.savefig("convergence_L16.pdf")
println("  已保存 convergence_L16.pdf")


# ==========================================
# 任务2：<m^2>随系统尺寸L的变化
# 对比 混合算法(Metropolis+SW) vs 纯Metropolis
# L = 4, 8, 16, 32
# ==========================================

function run_m2_vs_L(; L::Int, sigma::Float64, K::Float64=1.0,
                       therm_sweeps::Int=THERM_SWEEPS, measure_sweeps::Int=MEASURE_SWEEPS,
                       use_sw::Bool=true, step_size::Float64=1.0)
    h = zeros(Float64, L, L)

    # 热化
    for i in 1:therm_sweeps
        GeneralizedSOS.metropolis_sweep!(h, L, K, sigma, step_size)
        if use_sw && i % 5 == 0
            GeneralizedSOS.swendsen_wang_sweep!(h, L, K, sigma)
        end
    end

    # 测量
    roughnesses = Float64[]
    for i in 1:measure_sweeps
        GeneralizedSOS.metropolis_sweep!(h, L, K, sigma, step_size)
        if use_sw && i % 5 == 0
            GeneralizedSOS.swendsen_wang_sweep!(h, L, K, sigma)
        end
        push!(roughnesses, GeneralizedSOS.measure_roughness(h))
    end

    return mean(roughnesses), std(roughnesses) / sqrt(length(roughnesses))
end

println("\n=== 任务2：<m^2> vs L，对比混合算法与纯Metropolis ===")
L_list = [4, 8, 16, 32, 64, 128]

fig2, axes = subplots(1, 3, figsize=(14, 5))

for (idx, sigma) in enumerate(sigma_list)
    println("  sigma=$sigma ...")
    m2_hybrid  = Float64[]
    m2_metro   = Float64[]
    err_hybrid = Float64[]
    err_metro  = Float64[]

    for L in L_list
        println("    L=$L, hybrid ...")
        m, e = run_m2_vs_L(L=L, sigma=sigma, use_sw=true)
        push!(m2_hybrid, m); push!(err_hybrid, e)

        println("    L=$L, metro only ...")
        m, e = run_m2_vs_L(L=L, sigma=sigma, use_sw=false)
        push!(m2_metro, m); push!(err_metro, e)
    end

    ax = axes[idx]
    ax.errorbar(L_list, m2_hybrid, yerr=err_hybrid, marker="o", label="Metropolis + SW", linewidth=2, capsize=4)
    ax.errorbar(L_list, m2_metro,  yerr=err_metro,  marker="s", label="Metropolis only", linewidth=2, capsize=4, linestyle="--")
    ax.set_xlabel("System Size L", fontsize=12)
    ax.set_ylabel("\$\\langle m^2 \\rangle\$", fontsize=12)
    ax.set_title("\$\\sigma = $sigma\$", fontsize=13)
    ax.set_xscale("log", base=2)
    ax.set_xticks(L_list)
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.legend(fontsize=10)
    ax.grid(true, alpha=0.3)
end

fig2.suptitle("\$\\langle m^2 \\rangle\$ vs System Size \$L\$: Hybrid vs Pure Metropolis", fontsize=14)
fig2.tight_layout()
fig2.savefig("m2_vs_L_comparison.pdf")
println("  已保存 m2_vs_L_comparison.pdf")

println("\n全部完成。")
