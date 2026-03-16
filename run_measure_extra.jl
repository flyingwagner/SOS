include("SOSCore.jl")
include("Measure.jl")
using .SOSCore
using .Measure
using JLD2
using Statistics
using FFTW

# ==========================================
# 命令行参数解析
# ==========================================
function parse_args()
    args = ARGS
    if length(args) < 2
        println("用法: julia run_measure_extra.jl <L> <sigma> [选项]")
        println("  <L>      系统尺寸（整数）")
        println("  <sigma>  相互作用指数（浮点数）")
        println("选项:")
        println("  --therm    <N>   热化步数（默认 10000）")
        println("  --measure  <N>   测量步数（默认 100000）")
        println("  --binwidth <x>   histogram bin 宽度（默认 0.01）")
        println("  --maxdh    <x>   histogram 最大 |Δh|（默认 20.0）")
        println("  --no-sw          关闭 Swendsen-Wang 更新")
        println("  --datadir  <dir> 数据目录（默认 data_extra）")
        exit(1)
    end

    L     = parse(Int,     args[1])
    sigma = parse(Float64, args[2])

    therm     = 10000
    measure   = 100000
    bin_width = 0.01
    max_dh    = 20.0
    use_sw    = true
    data_dir  = "data_extra"

    i = 3
    while i <= length(args)
        if args[i] == "--therm"
            therm = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--measure"
            measure = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--binwidth"
            bin_width = parse(Float64, args[i+1]); i += 2
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

    return L, sigma, therm, measure, bin_width, max_dh, use_sw, data_dir
end

# ==========================================
# 主程序
# ==========================================
const K = 1.0

L, sigma, THERM, MEASURE, BIN_WIDTH, MAX_DH, USE_SW, DATA_DIR = parse_args()
NBINS = round(Int, MAX_DH / BIN_WIDTH)

mkpath(DATA_DIR)

key      = "L$(L)_s$(sigma)"
out_file = joinpath(DATA_DIR, "$(key)_extra.jld2")

println("========== L=$L, σ=$sigma ==========")
println("  热化: $THERM  测量: $MEASURE")
println("  histogram: binwidth=$BIN_WIDTH  maxdh=$MAX_DH  nbins=$NBINS")
println("  SW更新: $USE_SW  输出: $out_file")

h = zeros(Float64, L, L)

println("--- 热化中 ($THERM sweeps) ---")
thermalize!(h, L, K, sigma; therm_sweeps=THERM, use_sw=USE_SW)

# 累积 histogram
hist_counts = zeros(Float64, NBINS)

# 累积全 k 的 |h_k|²（对所有 k 求平均）
hk_power_sum = zeros(Float64, L, L)

println("--- 测量中 ($MEASURE sweeps) ---")
for i in 1:MEASURE
    metropolis_sweep!(h, L, K, sigma, 1.0)
    if USE_SW && i % 5 == 0
        swendsen_wang_sweep!(h, L, K, sigma)
    end

    # 累积 |h_i - h_j| histogram
    measure_height_diff_histogram!(hist_counts, h, L, BIN_WIDTH)

    # 累积 |h_k|²
    hk = measure_hk_full(h, L)
    hk_power_sum .+= abs2.(hk)

    if i % 10000 == 0
        println("  sweep $i / $MEASURE")
    end
end

# 归一化
hist_counts ./= MEASURE
hk_power_avg = hk_power_sum ./ MEASURE

# 提取 bin 中心
bin_centers = [(i - 0.5) * BIN_WIDTH for i in 1:NBINS]

jldsave(out_file;
    L, sigma, K,
    bin_centers, hist_counts, BIN_WIDTH,
    hk_power_avg,
)
println("\n========== 结果 ==========")
println("  histogram: $(sum(hist_counts)) 条键/sweep (理论值 $(2*L*L))")
println("  |h_k|² 平均值: $(mean(hk_power_avg))")
println("  已保存至 $out_file")
