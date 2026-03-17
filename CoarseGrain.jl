module CoarseGrain

export block_average_2x2!, define_shells, accumulate_kspace!, accumulate_realspace!

"""
    block_average_2x2!(h_out, h_in)

In-place 2×2 block averaging: M×M → M/2 × M/2.
"""
function block_average_2x2!(h_out::Matrix{Float64}, h_in::Matrix{Float64})
    M = size(h_in, 1)
    M2 = M ÷ 2
    @inbounds for j in 1:M2
        j2 = 2j - 1
        for i in 1:M2
            i2 = 2i - 1
            h_out[i, j] = (h_in[i2, j2] + h_in[i2+1, j2] +
                           h_in[i2, j2+1] + h_in[i2+1, j2+1]) * 0.25
        end
    end
end

"""
    define_shells(M; n_shells=20)

Define logarithmic radial k-shells for a M×M lattice.

Returns: (shell_id, shell_count, shell_k_mean, k_mag)
- shell_id[i,j] : shell index (0 for zero mode)
- shell_count[s] : number of k-points in shell s
- shell_k_mean[s]: mean |k| of shell s
- k_mag[i,j]    : |k| for each point
"""
function define_shells(M::Int; n_shells::Int=20)
    k_mag = zeros(Float64, M, M)
    for j in 1:M, i in 1:M
        nx = i <= M ÷ 2 + 1 ? i - 1 : i - 1 - M
        ny = j <= M ÷ 2 + 1 ? j - 1 : j - 1 - M
        k_mag[i, j] = sqrt((nx * 2π / M)^2 + (ny * 2π / M)^2)
    end

    k_nonzero = filter(>(0), vec(k_mag))
    k_min = minimum(k_nonzero)
    k_max = maximum(k_nonzero)

    log_min = log(k_min) - 0.01
    log_max = log(k_max) + 0.01
    log_step = (log_max - log_min) / n_shells

    shell_id = zeros(Int, M, M)
    for j in 1:M, i in 1:M
        km = k_mag[i, j]
        km == 0.0 && continue
        s = clamp(floor(Int, (log(km) - log_min) / log_step) + 1, 1, n_shells)
        shell_id[i, j] = s
    end

    shell_count  = zeros(Int, n_shells)
    shell_k_sum  = zeros(Float64, n_shells)
    for j in 1:M, i in 1:M
        s = shell_id[i, j]
        if s > 0
            shell_count[s] += 1
            shell_k_sum[s] += k_mag[i, j]
        end
    end
    shell_k_mean = [c > 0 ? shell_k_sum[s] / c : 0.0
                    for (s, c) in enumerate(shell_count)]

    return shell_id, shell_count, shell_k_mean, k_mag
end

"""
    accumulate_kspace!(hk2_acc, hk4_acc, S_sum, SS_sum,
                       hk, shell_id, shell_count, n_shells, S_buf)

Accumulate per-k-point |h_k|² and |h_k|⁴, plus shell-averaged power
S_α and cross-shell products S_α S_β, from one FFT snapshot `hk`.
"""
function accumulate_kspace!(hk2_acc::Matrix{Float64}, hk4_acc::Matrix{Float64},
                            S_sum::Vector{Float64}, SS_sum::Matrix{Float64},
                            hk::Matrix{ComplexF64},
                            shell_id::Matrix{Int}, shell_count::Vector{Int},
                            n_shells::Int, S_buf::Vector{Float64})
    M = size(hk, 1)
    fill!(S_buf, 0.0)

    @inbounds for j in 1:M, i in 1:M
        p2 = abs2(hk[i, j])
        hk2_acc[i, j] += p2
        hk4_acc[i, j] += p2 * p2
        s = shell_id[i, j]
        if s > 0
            S_buf[s] += p2
        end
    end

    # shell average: S_α = (1/N_α) Σ_{k∈α} |h_k|²
    @inbounds for s in 1:n_shells
        if shell_count[s] > 0
            S_buf[s] /= shell_count[s]
        end
    end

    @inbounds begin
        for s in 1:n_shells
            S_sum[s] += S_buf[s]
        end
        for β in 1:n_shells, α in 1:n_shells
            SS_sum[α, β] += S_buf[α] * S_buf[β]
        end
    end
end

"""
    accumulate_realspace!(hist, moms, h, M, bin_width, nbins)

Accumulate nearest-neighbor |Δh̄| histogram and moments
[Σ(Δh)², Σ(Δh)⁴, Σ(Δh)⁶] on an M×M lattice with PBC.
Each bond (right + up) counted once.
"""
function accumulate_realspace!(hist::Vector{Float64}, moms::Vector{Float64},
                               h::Matrix{Float64}, M::Int,
                               bin_width::Float64, nbins::Int)
    @inbounds for j in 1:M, i in 1:M
        # right neighbor
        ni = i == M ? 1 : i + 1
        dh  = h[i, j] - h[ni, j]
        adh = abs(dh)
        bin = clamp(floor(Int, adh / bin_width) + 1, 1, nbins)
        hist[bin] += 1.0
        dh2 = dh * dh
        moms[1] += dh2
        moms[2] += dh2 * dh2
        moms[3] += dh2 * dh2 * dh2

        # up neighbor
        nj = j == M ? 1 : j + 1
        dh  = h[i, j] - h[i, nj]
        adh = abs(dh)
        bin = clamp(floor(Int, adh / bin_width) + 1, 1, nbins)
        hist[bin] += 1.0
        dh2 = dh * dh
        moms[1] += dh2
        moms[2] += dh2 * dh2
        moms[3] += dh2 * dh2 * dh2
    end
end

end # module
