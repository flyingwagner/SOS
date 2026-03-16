module Measure

using Statistics
using FFTW

export measure_roughness, measure_structure_factor_kmin,
       measure_height_diff_histogram!, measure_hk_full

# ==========================================
# 现有测量量
# ==========================================

"""
    measure_roughness(h) -> Float64

W² = (1/N) Σ h_i²（零模已消除）
"""
function measure_roughness(h::Matrix{Float64})
    return mean(h .^ 2)
end

"""
    measure_structure_factor_kmin(h, L) -> (|h_k|², Re(h_k), Im(h_k))

k_min = (2π/L, 0) 处的归一化傅里叶模式。
返回 (|h_k|²/N, Re(h_k)/√N, Im(h_k)/√N)。
"""
function measure_structure_factor_kmin(h::Matrix{Float64}, L::Int)
    sum_cos = 0.0
    sum_sin = 0.0
    k_x = 2.0 * π / L
    for y in 1:L
        for x in 1:L
            phase = k_x * x
            sum_cos += h[x, y] * cos(phase)
            sum_sin += h[x, y] * sin(phase)
        end
    end
    sqrtN = sqrt(L * L)
    re_hk = sum_cos / sqrtN
    im_hk = sum_sin / sqrtN
    return re_hk^2 + im_hk^2, re_hk, im_hk
end

# ==========================================
# 新增测量量
# ==========================================

"""
    measure_height_diff_histogram!(counts, h, L, bin_width) -> nothing

累积 |h_i - h_j| 的 histogram 到 counts 向量中（in-place 累加）。
每条最近邻键只计一次（只遍历右邻居和上邻居）。
bin 编号：floor(|Δh| / bin_width) + 1，超出范围的归入最后一个 bin。
"""
function measure_height_diff_histogram!(counts::Vector{Float64}, h::Matrix{Float64},
                                        L::Int, bin_width::Float64)
    nbins = length(counts)
    for y in 1:L
        for x in 1:L
            # 右邻居
            nx = x == L ? 1 : x + 1
            dh = abs(h[x, y] - h[nx, y])
            bin = clamp(floor(Int, dh / bin_width) + 1, 1, nbins)
            counts[bin] += 1.0

            # 上邻居
            ny = y == L ? 1 : y + 1
            dh = abs(h[x, y] - h[x, ny])
            bin = clamp(floor(Int, dh / bin_width) + 1, 1, nbins)
            counts[bin] += 1.0
        end
    end
end

"""
    measure_hk_full(h, L) -> Matrix{ComplexF64}

用 FFT 计算全部 k 的 h_k = FFT(h) / √N。
返回 L×L 的复数矩阵。
"""
function measure_hk_full(h::Matrix{Float64}, L::Int)
    return fft(h) / sqrt(L * L)
end

end # module
