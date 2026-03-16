include("SOSCore.jl")
include("Measure.jl")
include("BlockStats.jl")
using .SOSCore
using .Measure
using .BlockStats
using JLD2
using Statistics

# ==========================================
# 模拟参数
# ==========================================
const L_LIST      = [16]
const SIGMA_LIST  = [0.5, 1.0, 2.0, 4.0]
const K           = 1.0
const THERM       = 10000
const MEASURE     = 100000
const USE_SW      = true
const DATA_DIR    = "data"
const BLOCK_SIZE  = 100

mkpath(DATA_DIR)

# ==========================================
# 主循环：遍历所有 (L, sigma) 组合
# ==========================================
for sigma in SIGMA_LIST
    for L in L_LIST
        key = "L$(L)_s$(sigma)"
        out_file = joinpath(DATA_DIR, "$(key).jld2")
        println("\n========== L=$L, σ=$sigma ==========")

        h = zeros(Float64, L, L)

        println("--- 热化中 ---")
        thermalize!(h, L, K, sigma; therm_sweeps=THERM, use_sw=USE_SW)

        W2_series    = Vector{Float64}(undef, MEASURE)
        chi_series   = Vector{Float64}(undef, MEASURE)
        rehk_series  = Vector{Float64}(undef, MEASURE)
        imhk_series  = Vector{Float64}(undef, MEASURE)

        println("--- 测量中 ---")
        for i in 1:MEASURE
            metropolis_sweep!(h, L, K, sigma, 1.0)
            if USE_SW && i % 5 == 0
                swendsen_wang_sweep!(h, L, K, sigma)
            end
            W2_series[i] = measure_roughness(h)
            hk2, re, im  = measure_structure_factor_kmin(h, L)
            chi_series[i]  = hk2
            rehk_series[i] = re
            imhk_series[i] = im
        end

        W2_mean,  W2_err,  W2_rho1,  W2_bs  = block_stats(W2_series,  BLOCK_SIZE)
        chi_mean, chi_err, chi_rho1, chi_bs = block_stats(chi_series, BLOCK_SIZE)

        # chi' = <|h_k|^2> - |<h_k>|^2
        n_blocks = div(MEASURE, BLOCK_SIZE)
        chip_blocks = Vector{Float64}(undef, n_blocks)
        for b in 1:n_blocks
            rng = (b-1)*BLOCK_SIZE+1 : b*BLOCK_SIZE
            chip_blocks[b] = mean(chi_series[rng]) - mean(rehk_series[rng])^2 - mean(imhk_series[rng])^2
        end
        chip_mean, chip_err, chip_rho1, chip_bs = block_stats(chip_blocks, 1)

        jldsave(out_file;
            L, sigma, K,
            W2_mean,   W2_err,   W2_rho1,   W2_bs,
            chi_mean,  chi_err,  chi_rho1,  chi_bs,
            chip_mean, chip_err, chip_rho1, chip_bs,
            W2_series, chi_series, rehk_series, imhk_series,
        )
        println("  W²  = $(round(W2_mean,  digits=4)) ± $(round(W2_err,  digits=4))  ρ(1)=$(round(W2_rho1,  digits=4))  block=$(W2_bs)")
        println("  χ   = $(round(chi_mean,  digits=4)) ± $(round(chi_err,  digits=4))  ρ(1)=$(round(chi_rho1,  digits=4))  block=$(chi_bs)")
        println("  χ'  = $(round(chip_mean, digits=4)) ± $(round(chip_err, digits=4))  ρ(1)=$(round(chip_rho1, digits=4))  block=$(chip_bs*BLOCK_SIZE)")
        println("  已保存至 $out_file")
    end
end
