include("SOSCore3D.jl")
include("Measure3D.jl")
include("BlockStats.jl")
using .SOSCore3D
using .Measure3D
using .BlockStats
using JLD2
using Statistics
using FFTW

# ==========================================
# 3D SOS 主驱动：模拟 3 维立方格子上的 generalized SOS。
# H = K Σ_{<ij>} |h_i - h_j|^σ，零模约束 Σ h = 0，PBC，6 个最近邻。
#
# 主要观测量 M² = (1/N) Σ h_i²（与 report 中 notation 一致）。
# 内部变量名仍为 W2_* 以与 2D run.jl 兼容（旧数据/分析脚本字段名）。
# ==========================================

function parse_args()
    args = ARGS
    if length(args) < 2
        println("用法: julia run3d.jl <L> <sigma> [选项]")
        println("  <L>      系统线性尺寸（N = L³）")
        println("  <sigma>  相互作用指数")
        println("选项:")
        println("  --therm   <N>   热化步数（默认 10000）")
        println("  --measure <N>   测量步数（默认 100000）")
        println("  --block   <N>   blocking 初始块大小（默认 100）")
        println("  --nbins   <N>   histogram bin 数（默认 200）")
        println("  --maxdh   <x>   histogram 最大 |Δh|（默认 10.0）")
        println("  --no-sw         关闭 Swendsen-Wang 更新")
        println("  --datadir <dir> 数据目录（默认 data3d）")
        exit(1)
    end

    L     = parse(Int,     args[1])
    sigma = parse(Float64, args[2])

    therm    = 10000
    measure  = 100000
    block_sz = 100
    nbins    = 200
    max_dh   = 10.0
    use_sw   = true
    data_dir = "data3d"

    i = 3
    while i <= length(args)
        if args[i] == "--therm"
            therm = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--measure"
            measure = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--block"
            block_sz = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--nbins"
            nbins = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--maxdh"
            max_dh = parse(Float64, args[i+1]); i += 2
        elseif args[i] == "--no-sw"
            use_sw = false; i += 1
        elseif args[i] == "--datadir"
            data_dir = args[i+1]; i += 2
        else
            println("未知参数: $(args[i])"); exit(1)
        end
    end

    return L, sigma, therm, measure, block_sz, nbins, max_dh, use_sw, data_dir
end

const K = 1.0

L, sigma, THERM, MEASURE, BLOCK_SIZE, NBINS, MAX_DH, USE_SW, DATA_DIR = parse_args()
BIN_WIDTH = MAX_DH / NBINS

mkpath(DATA_DIR)

key      = "L$(L)_s$(sigma)_3d"
out_file = joinpath(DATA_DIR, "$(key).jld2")

println("========== 3D L=$L, σ=$sigma  (N=$(L^3)) ==========")
println("  热化: $THERM  测量: $MEASURE  block: $BLOCK_SIZE")
println("  histogram: nbins=$NBINS  max_dh=$MAX_DH  bin_width=$(round(BIN_WIDTH, digits=4))")
println("  SW更新: $USE_SW  输出: $out_file")

# ---- 初始化 ----
h = zeros(Float64, L, L, L)

if isfile(out_file)
    println("  检测到已有数据文件，尝试热启动...")
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

# ---- 热化 ----
if THERM > 0
    println("--- 热化中 ($THERM sweeps) ---")
    thermalize!(h, L, K, sigma; therm_sweeps=THERM, use_sw=USE_SW)
end

# ---- 测量缓冲区 ----
W2_series   = Vector{Float64}(undef, MEASURE)
chi_series  = Vector{Float64}(undef, MEASURE)
rehk_series = Vector{Float64}(undef, MEASURE)
imhk_series = Vector{Float64}(undef, MEASURE)

hist_counts    = zeros(Float64, NBINS)
hk_power_sum   = zeros(Float64, L, L, L)

# ---- 测量循环 ----
println("--- 测量中 ($MEASURE sweeps) ---")
for i in 1:MEASURE
    metropolis_sweep!(h, L, K, sigma, 1.0)
    if USE_SW && i % 5 == 0
        swendsen_wang_sweep!(h, L, K, sigma)
    end

    W2_series[i]   = measure_roughness(h)
    hk2, re, im    = measure_structure_factor_kmin(h, L)
    chi_series[i]  = hk2
    rehk_series[i] = re
    imhk_series[i] = im

    measure_height_diff_histogram!(hist_counts, h, L, BIN_WIDTH)

    hk = measure_hk_full(h, L)
    hk_power_sum .+= abs2.(hk)

    if i % 10000 == 0
        println("  sweep $i / $MEASURE  M²=$(round(W2_series[i], digits=4))")
    end
end

# ---- 统计分析 ----
W2_mean,  W2_err,  W2_rho1,  W2_bs  = block_stats(W2_series,  BLOCK_SIZE)
chi_mean, chi_err, chi_rho1, chi_bs = block_stats(chi_series, BLOCK_SIZE)

n_blocks = div(MEASURE, BLOCK_SIZE)
chip_blocks = Vector{Float64}(undef, n_blocks)
for b in 1:n_blocks
    rng = (b-1)*BLOCK_SIZE+1 : b*BLOCK_SIZE
    chip_blocks[b] = mean(chi_series[rng]) - mean(rehk_series[rng])^2 - mean(imhk_series[rng])^2
end
chip_mean, chip_err, chip_rho1, chip_bs = block_stats(chip_blocks, 1)

hist_counts ./= MEASURE
bin_centers = [(i - 0.5) * BIN_WIDTH for i in 1:NBINS]

hk_power_avg = hk_power_sum ./ MEASURE

# ---- 保存 ----
jldsave(out_file;
    L, sigma, K, dim = 3,
    W2_mean,   W2_err,   W2_rho1,   W2_bs,
    chi_mean,  chi_err,  chi_rho1,  chi_bs,
    chip_mean, chip_err, chip_rho1, chip_bs,
    W2_series, chi_series, rehk_series, imhk_series,
    bin_centers, hist_counts, BIN_WIDTH,
    hk_power_avg,
    h_last = copy(h),
)

println("\n========== 结果 ==========")
println("  M²   = $(round(W2_mean,  digits=5)) ± $(round(W2_err,  digits=5))  ρ(1)=$(round(W2_rho1,  digits=4))  block=$(W2_bs)")
println("  χ    = $(round(chi_mean,  digits=5)) ± $(round(chi_err,  digits=5))  ρ(1)=$(round(chi_rho1,  digits=4))  block=$(chi_bs)")
println("  χ'   = $(round(chip_mean, digits=5)) ± $(round(chip_err, digits=5))  ρ(1)=$(round(chip_rho1, digits=4))  block=$(chip_bs*BLOCK_SIZE)")
println("  histogram: $(sum(hist_counts)) 键/sweep（理论 $(3*L*L*L)）")
println("  |h_k|² 全谱均值: $(round(mean(hk_power_avg), digits=5))")
println("  已保存至 $out_file")
