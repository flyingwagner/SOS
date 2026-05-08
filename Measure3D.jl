module Measure3D

using Statistics
using FFTW

export measure_roughness, measure_correlation_x,
       measure_structure_factor_kmin, measure_height_diff_histogram!,
       measure_hk_full

"""
    measure_roughness(h) -> Float64

M² = (1/N) Σ h_i² ，零模已通过 zero_mode_shift! 消除。
"""
function measure_roughness(h::Array{Float64,3})
    return mean(h .^ 2)
end

"""
    measure_correlation_x(h, L) -> Vector{Float64}

沿 x 方向的高度差平方关联 G(r) = <(h(x+r,y,z) - h(x,y,z))²>，r = 1..L/2。
对所有 (y,z) 平面平均。
"""
function measure_correlation_x(h::Array{Float64,3}, L::Int)
    max_r = div(L, 2)
    G = zeros(Float64, max_r)
    invN = 1.0 / (L * L * L)
    for r in 1:max_r
        s = 0.0
        for z in 1:L, y in 1:L, x in 1:L
            nx = mod1(x + r, L)
            s += (h[nx, y, z] - h[x, y, z])^2
        end
        G[r] = s * invN
    end
    return G
end

"""
    measure_structure_factor_kmin(h, L) -> (|h_k|², Re(h_k), Im(h_k))

k_min = (2π/L, 0, 0) 处归一化 Fourier 模。返回 (|h_k|²/N, Re(h_k)/√N, Im(h_k)/√N)。
"""
function measure_structure_factor_kmin(h::Array{Float64,3}, L::Int)
    sum_cos = 0.0
    sum_sin = 0.0
    k_x = 2.0 * π / L
    for z in 1:L, y in 1:L, x in 1:L
        phase = k_x * x
        sum_cos += h[x, y, z] * cos(phase)
        sum_sin += h[x, y, z] * sin(phase)
    end
    sqrtN = sqrt(L * L * L)
    re_hk = sum_cos / sqrtN
    im_hk = sum_sin / sqrtN
    return re_hk^2 + im_hk^2, re_hk, im_hk
end

"""
    measure_height_diff_histogram!(counts, h, L, bin_width) -> nothing

|h_i - h_j| 的 histogram。每条最近邻键只计一次：遍历 +x, +y, +z 三个方向。
"""
function measure_height_diff_histogram!(counts::Vector{Float64}, h::Array{Float64,3},
                                        L::Int, bin_width::Float64)
    nbins = length(counts)
    for z in 1:L, y in 1:L, x in 1:L
        hi = h[x, y, z]

        nx = x == L ? 1 : x + 1
        dh = abs(hi - h[nx, y, z])
        bin = clamp(floor(Int, dh / bin_width) + 1, 1, nbins)
        counts[bin] += 1.0

        ny = y == L ? 1 : y + 1
        dh = abs(hi - h[x, ny, z])
        bin = clamp(floor(Int, dh / bin_width) + 1, 1, nbins)
        counts[bin] += 1.0

        nz = z == L ? 1 : z + 1
        dh = abs(hi - h[x, y, nz])
        bin = clamp(floor(Int, dh / bin_width) + 1, 1, nbins)
        counts[bin] += 1.0
    end
end

"""
    measure_hk_full(h, L) -> Array{ComplexF64,3}

3D FFT，h_k = FFT(h) / √N，N = L³。
"""
function measure_hk_full(h::Array{Float64,3}, L::Int)
    return fft(h) / sqrt(L * L * L)
end

end # module
