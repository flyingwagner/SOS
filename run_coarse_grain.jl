include("SOSCore.jl")
include("CoarseGrain.jl")
using .SOSCore
using .CoarseGrain
using JLD2
using Statistics
using FFTW
using LinearAlgebra: mul!

# ==========================================
# 命令行参数
# ==========================================
function parse_args()
    args = ARGS
    if length(args) < 2
        println("用法: julia run_coarse_grain.jl <L> <sigma> [选项]")
        println("  <L>      系统尺寸（2 的幂次）")
        println("  <sigma>  相互作用指数")
        println("选项:")
        println("  --therm    <N>   热化步数（默认 10000）")
        println("  --measure  <N>   测量步数（默认 100000）")
        println("  --no-sw          关闭 Swendsen-Wang 更新")
        println("  --datadir  <dir> 数据目录（默认 data_cg）")
        exit(1)
    end
    L     = parse(Int, args[1])
    sigma = parse(Float64, args[2])
    therm = 10000; measure = 100000; use_sw = true; data_dir = "data_cg"
    i = 3
    while i <= length(args)
        if     args[i] == "--therm";   therm    = parse(Int, args[i+1]);     i += 2
        elseif args[i] == "--measure"; measure  = parse(Int, args[i+1]);     i += 2
        elseif args[i] == "--no-sw";   use_sw   = false;                     i += 1
        elseif args[i] == "--datadir"; data_dir = args[i+1];                 i += 2
        else   println("未知参数: $(args[i])"); exit(1)
        end
    end
    return L, sigma, therm, measure, use_sw, data_dir
end

# ==========================================
# 主程序
# ==========================================
const K = 1.0
L, sigma, THERM, MEASURE, USE_SW, DATA_DIR = parse_args()
@assert ispow2(L) "L must be a power of 2"

# CG 层级: b = 1, 2, 4, ... 自适应确定最大 block
# 实空间最小 CG 格子 M_min=8，k 空间最小 M_min=16
const MIN_M_REAL = 8
const MIN_M_K    = 16
const MAX_B_REAL = L ÷ MIN_M_REAL       # 实空间最大 block
const MAX_B_K    = L ÷ MIN_M_K          # k 空间最大 block

# 全部层（实空间）: b = 1, 2, 4, ..., MAX_B_REAL
const BLOCK_SIZES = [1 << k for k in 0:Int(log2(MAX_B_REAL))]
const N_LEVELS    = length(BLOCK_SIZES)
const CG_M        = [L ÷ b for b in BLOCK_SIZES]

# k 空间层: b ≤ MAX_B_K
const N_LEVELS_K  = count(b -> b <= MAX_B_K, BLOCK_SIZES)
const N_SHELLS    = 20

# Histogram 参数
const BIN_WIDTH = 0.01
const MAX_DH    = 20.0
const NBINS     = round(Int, MAX_DH / BIN_WIDTH)

mkpath(DATA_DIR)
out_file = joinpath(DATA_DIR, "L$(L)_s$(sigma)_cg.jld2")

println("========== L=$L, σ=$sigma ==========")
println("  热化=$THERM  测量=$MEASURE  SW=$USE_SW")
println("  CG 层级: b=$(BLOCK_SIZES), M=$(CG_M)")
println("  输出: $out_file")

# ==========================================
# 预计算 shell 定义（每个 k 空间层）
# ==========================================
shell_ids    = Vector{Matrix{Int}}(undef, N_LEVELS_K)
shell_counts = Vector{Vector{Int}}(undef, N_LEVELS_K)
shell_k_means = Vector{Vector{Float64}}(undef, N_LEVELS_K)
k_mags       = Vector{Matrix{Float64}}(undef, N_LEVELS_K)

for lev in 1:N_LEVELS_K
    M = CG_M[lev]
    shell_ids[lev], shell_counts[lev], shell_k_means[lev], k_mags[lev] =
        define_shells(M; n_shells=N_SHELLS)
    n_active = count(>(0), shell_counts[lev])
    println("  level $lev (b=$(BLOCK_SIZES[lev]), M=$M): $n_active active shells")
end

# ==========================================
# 预分配累积量
# ==========================================
# --- k 空间 per-k-point ---
hk2_sum = [zeros(Float64, CG_M[lev], CG_M[lev]) for lev in 1:N_LEVELS_K]
hk4_sum = [zeros(Float64, CG_M[lev], CG_M[lev]) for lev in 1:N_LEVELS_K]

# --- k 空间 shell 关联 ---
shell_S_sum  = [zeros(Float64, N_SHELLS) for _ in 1:N_LEVELS_K]
shell_SS_sum = [zeros(Float64, N_SHELLS, N_SHELLS) for _ in 1:N_LEVELS_K]

# --- 实空间 histogram + moments ---
hist_all = [zeros(Float64, NBINS) for _ in 1:N_LEVELS]
moms_all = [zeros(Float64, 3) for _ in 1:N_LEVELS]   # [Σ(Δh)², Σ(Δh)⁴, Σ(Δh)⁶]

# --- W² 和 χ_kmin（交叉验证）---
W2_sum   = Ref(0.0)
chi_sum  = Ref(0.0)

# ==========================================
# 工作缓冲区
# ==========================================
# CG 场缓冲区（层级式复用）
h_cg = [Matrix{Float64}(undef, CG_M[lev], CG_M[lev]) for lev in 2:N_LEVELS]

# FFT 缓冲区（b=1 用预分配 + plan，CG 层用 fft()）
const hk_in  = Matrix{ComplexF64}(undef, L, L)
const hk_out = Matrix{ComplexF64}(undef, L, L)
const fft_plan = plan_fft(hk_in)
const sqrtN = sqrt(Float64(L * L))

# Shell 工作缓冲区
S_buf = zeros(Float64, N_SHELLS)

# ==========================================
# 高度场初始化 + 热启动
# ==========================================
h = zeros(Float64, L, L)

if isfile(out_file)
    try
        jldopen(out_file, "r") do f
            if haskey(f, "h_last")
                h .= f["h_last"]
                println("  热启动成功，跳过热化。")
                global THERM = 0
            end
        end
    catch e
        println("  热启动失败（$e），从零开始。")
    end
end

# ==========================================
# 热化
# ==========================================
if THERM > 0
    println("--- 热化中 ($THERM sweeps) ---")
    thermalize!(h, L, K, sigma; therm_sweeps=THERM, use_sw=USE_SW)
end

# ==========================================
# 测量循环
# ==========================================
println("--- 测量中 ($MEASURE sweeps) ---")
t_start = time()

for sweep in 1:MEASURE
    # ---- MC 更新 ----
    metropolis_sweep!(h, L, K, sigma, 1.0)
    if USE_SW && sweep % 5 == 0
        swendsen_wang_sweep!(h, L, K, sigma)
    end

    # ---- Level 1 (b=1): 原始场 ----
    # 实空间
    accumulate_realspace!(hist_all[1], moms_all[1], h, L, BIN_WIDTH, NBINS)

    # W²
    w2 = 0.0
    @inbounds for j in 1:L, i in 1:L
        w2 += h[i, j]^2
    end
    W2_sum[] += w2 / (L * L)

    # k 空间（预分配 FFT）
    hk_in .= h
    mul!(hk_out, fft_plan, hk_in)
    hk_out ./= sqrtN

    # χ_kmin = |h_{kmin}|²，kmin = (2π/L, 0) → index (2, 1)
    chi_sum[] += abs2(hk_out[2, 1])

    accumulate_kspace!(hk2_sum[1], hk4_sum[1],
                       shell_S_sum[1], shell_SS_sum[1],
                       hk_out, shell_ids[1], shell_counts[1],
                       N_SHELLS, S_buf)

    # ---- 层级 CG: level 2–8 ----
    h_prev = h
    for lev in 2:N_LEVELS
        block_average_2x2!(h_cg[lev - 1], h_prev)
        hf = h_cg[lev - 1]
        M  = CG_M[lev]

        # 实空间
        accumulate_realspace!(hist_all[lev], moms_all[lev], hf, M, BIN_WIDTH, NBINS)

        # k 空间（仅 lev ≤ 7，即 b ≤ 64，M ≥ 16）
        if lev <= N_LEVELS_K
            sqrtM2 = sqrt(Float64(M * M))
            hk_cg  = fft(hf)
            hk_cg ./= sqrtM2
            accumulate_kspace!(hk2_sum[lev], hk4_sum[lev],
                               shell_S_sum[lev], shell_SS_sum[lev],
                               hk_cg, shell_ids[lev], shell_counts[lev],
                               N_SHELLS, S_buf)
        end

        h_prev = hf
    end

    # ---- 进度报告 ----
    if sweep % 10000 == 0
        elapsed = time() - t_start
        rate = sweep / elapsed
        eta  = (MEASURE - sweep) / rate
        println("  sweep $sweep / $MEASURE  " *
                "($(round(rate, digits=1)) sweeps/s, " *
                "ETA $(round(eta / 60, digits=1)) min)")
    end
end

elapsed_total = time() - t_start
println("--- 测量完成 ($(round(elapsed_total / 60, digits=1)) min) ---")

# ==========================================
# 归一化
# ==========================================
for lev in 1:N_LEVELS_K
    hk2_sum[lev] ./= MEASURE
    hk4_sum[lev] ./= MEASURE
    shell_S_sum[lev]  ./= MEASURE
    shell_SS_sum[lev] ./= MEASURE
end

for lev in 1:N_LEVELS
    hist_all[lev] ./= MEASURE
end

bin_centers     = [(i - 0.5) * BIN_WIDTH for i in 1:NBINS]
bonds_per_sweep = [2 * CG_M[lev]^2 for lev in 1:N_LEVELS]
total_bonds     = [MEASURE * bps for bps in bonds_per_sweep]

# 矩归一化 → 均值
mean_dh2 = [moms_all[lev][1] / total_bonds[lev] for lev in 1:N_LEVELS]
mean_dh4 = [moms_all[lev][2] / total_bonds[lev] for lev in 1:N_LEVELS]
mean_dh6 = [moms_all[lev][3] / total_bonds[lev] for lev in 1:N_LEVELS]
kurtosis = [mean_dh4[lev] / mean_dh2[lev]^2 - 3.0 for lev in 1:N_LEVELS]

W2_avg  = W2_sum[]  / MEASURE
chi_avg = chi_sum[] / MEASURE

# ==========================================
# 输出摘要
# ==========================================
println("\n========== 结果摘要 ==========")
println("  W²  = $(round(W2_avg, digits=6))")
println("  χ_k = $(round(chi_avg, digits=6))")
println()
println("  b       M     ⟨(Δh)²⟩      κ₄(excess)")
println("  " * "-"^50)
for lev in 1:N_LEVELS
    b = BLOCK_SIZES[lev]
    M = CG_M[lev]
    println("  $(lpad(b, 5))  $(lpad(M, 5))   " *
            "$(lpad(round(mean_dh2[lev], digits=6), 10))   " *
            "$(lpad(round(kurtosis[lev], digits=4), 8))")
end

# ==========================================
# 保存
# ==========================================
jldsave(out_file;
    # 元数据
    L, sigma, K, MEASURE, BIN_WIDTH, MAX_DH, NBINS, N_SHELLS,
    block_sizes     = collect(BLOCK_SIZES),
    cg_M            = collect(CG_M),
    elapsed_seconds = elapsed_total,
    # k 空间 per-k-point（已除以 MEASURE）
    hk2_avg   = hk2_sum,        # Vector{Matrix{Float64}}, length N_LEVELS_K
    hk4_avg   = hk4_sum,
    k_mags    = k_mags,         # |k| for each (i,j)
    # k 空间 shell 关联（已除以 MEASURE）
    shell_S_avg     = shell_S_sum,
    shell_SS_avg    = shell_SS_sum,
    shell_k_means   = shell_k_means,
    shell_counts    = shell_counts,
    # 实空间
    hist            = hist_all,          # Vector{Vector{Float64}}, length N_LEVELS
    bin_centers     = bin_centers,
    moms_raw        = moms_all,          # 原始累积矩（未归一化）
    bonds_per_sweep = bonds_per_sweep,
    mean_dh2        = mean_dh2,
    mean_dh4        = mean_dh4,
    mean_dh6        = mean_dh6,
    kurtosis_excess = kurtosis,
    # 交叉验证
    W2_avg  = W2_avg,
    chi_avg = chi_avg,
    # 最终构型（供热启动）
    h_last = copy(h),
)
println("\n  已保存至 $out_file")
