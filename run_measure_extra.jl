include("SOSCore.jl")
include("Measure.jl")
using .SOSCore
using .Measure
using JLD2
using Statistics
using FFTW

# ==========================================
# 测量参数
# ==========================================
const L_LIST      = [64]
const SIGMA_LIST  = [0.5, 1.0, 2.0, 4.0]
const K           = 1.0
const THERM       = 10000
const MEASURE     = 100000   # histogram 收敛较快，不需要太多 sweep
const USE_SW      = true
const DATA_DIR    = "data_extra"

# Histogram 参数
const BIN_WIDTH   = 0.01
const MAX_DH      = 20.0
const NBINS       = round(Int, MAX_DH / BIN_WIDTH)

mkpath(DATA_DIR)

# ==========================================
# 主循环
# ==========================================
for sigma in SIGMA_LIST
    for L in L_LIST
        key = "L$(L)_s$(sigma)"
        out_file = joinpath(DATA_DIR, "$(key)_extra.jld2")
        println("\n========== L=$L, σ=$sigma ==========")

        h = zeros(Float64, L, L)

        println("--- 热化中 ---")
        thermalize!(h, L, K, sigma; therm_sweeps=THERM, use_sw=USE_SW)

        # 累积 histogram
        hist_counts = zeros(Float64, NBINS)

        # 累积全 k 的 |h_k|²（对所有 k 求平均）
        hk_power_sum = zeros(Float64, L, L)

        println("--- 测量中 ---")
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
        println("  已保存至 $out_file")
        println("  histogram: $(sum(hist_counts)) 条键/sweep (理论值 $(2*L*L))")
        println("  |h_k|² 平均值: $(mean(hk_power_avg))")
    end
end

println("\n全部完成。")
