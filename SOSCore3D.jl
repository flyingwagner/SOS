module SOSCore3D

using Random
using Statistics

export zero_mode_shift!, get_neighbors, metropolis_sweep!, swendsen_wang_sweep!, thermalize!

# ==========================================
# 并查集 (Union-Find)
# ==========================================
mutable struct UnionFind
    parent::Vector{Int}
    UnionFind(n::Int) = new(collect(1:n))
end

function find!(uf::UnionFind, i::Int)
    if uf.parent[i] == i
        return i
    end
    uf.parent[i] = find!(uf, uf.parent[i])
    return uf.parent[i]
end

function union!(uf::UnionFind, i::Int, j::Int)
    root_i = find!(uf, i)
    root_j = find!(uf, j)
    if root_i != root_j
        uf.parent[root_i] = root_j
    end
end

# ==========================================
# 基础工具：3D 立方格子，周期边界
# ==========================================
function zero_mode_shift!(h::Array{Float64,3})
    h .-= mean(h)
end

@inline function get_neighbors(x::Int, y::Int, z::Int, L::Int)
    xp = x == L ? 1 : x + 1
    xm = x == 1 ? L : x - 1
    yp = y == L ? 1 : y + 1
    ym = y == 1 ? L : y - 1
    zp = z == L ? 1 : z + 1
    zm = z == 1 ? L : z - 1
    return (xp, y, z), (xm, y, z), (x, yp, z), (x, ym, z), (x, y, zp), (x, y, zm)
end

@inline lin_idx(x::Int, y::Int, z::Int, L::Int) = ((z - 1) * L + (y - 1)) * L + x

# ==========================================
# Metropolis 单步
# ==========================================
function metropolis_sweep!(h::Array{Float64,3}, L::Int, K::Float64, sigma::Float64, step_size::Float64)
    accepted = 0
    N = L * L * L
    for _ in 1:N
        x, y, z = rand(1:L), rand(1:L), rand(1:L)
        old_h = h[x, y, z]
        new_h = old_h + step_size * (2.0 * rand() - 1.0)

        delta_E = 0.0
        for (nx, ny, nz) in get_neighbors(x, y, z, L)
            neighbor_h = h[nx, ny, nz]
            delta_E += K * (abs(new_h - neighbor_h)^sigma - abs(old_h - neighbor_h)^sigma)
        end

        if delta_E <= 0.0 || rand() < exp(-delta_E)
            h[x, y, z] = new_h
            accepted += 1
        end
    end
    zero_mode_shift!(h)
    return accepted / N
end

# ==========================================
# Swendsen-Wang 嵌入反射集团（3D）
# 在 +x, +y, +z 三个方向上枚举键，避免重复计数
# ==========================================
function swendsen_wang_sweep!(h::Array{Float64,3}, L::Int, K::Float64, sigma::Float64)
    N = L * L * L
    uf = UnionFind(N)
    h_ref = h[rand(1:L), rand(1:L), rand(1:L)]

    @inline function try_bond!(uf, i_idx, j_idx, hi, hj)
        delta_E = K * abs(2.0 * h_ref - hi - hj)^sigma - K * abs(hi - hj)^sigma
        if delta_E > 0.0 && rand() < 1.0 - exp(-delta_E)
            union!(uf, i_idx, j_idx)
        end
    end

    for z in 1:L, y in 1:L, x in 1:L
        i_idx = lin_idx(x, y, z, L)
        hi = h[x, y, z]
        si = hi > h_ref ? 1 : -1

        # +x
        nx = x == L ? 1 : x + 1
        hj = h[nx, y, z]
        sj = hj > h_ref ? 1 : -1
        if si == sj
            try_bond!(uf, i_idx, lin_idx(nx, y, z, L), hi, hj)
        end

        # +y
        ny = y == L ? 1 : y + 1
        hj = h[x, ny, z]
        sj = hj > h_ref ? 1 : -1
        if si == sj
            try_bond!(uf, i_idx, lin_idx(x, ny, z, L), hi, hj)
        end

        # +z
        nz = z == L ? 1 : z + 1
        hj = h[x, y, nz]
        sj = hj > h_ref ? 1 : -1
        if si == sj
            try_bond!(uf, i_idx, lin_idx(x, y, nz, L), hi, hj)
        end
    end

    cluster_flip = Dict{Int, Bool}()
    for i in 1:N
        root = find!(uf, i)
        if !haskey(cluster_flip, root)
            cluster_flip[root] = rand() < 0.5
        end
        if cluster_flip[root]
            # 解码 i -> (x, y, z)
            i0 = i - 1
            cx = mod(i0, L) + 1
            cy = mod(div(i0, L), L) + 1
            cz = div(i0, L * L) + 1
            h[cx, cy, cz] = 2.0 * h_ref - h[cx, cy, cz]
        end
    end
    zero_mode_shift!(h)
end

# ==========================================
# 热化
# ==========================================
"""
    thermalize!(h, L, K, sigma; therm_sweeps, use_sw, step_size)

3D 热化：therm_sweeps 步 Metropolis，每 5 步插入一次 SW。
"""
function thermalize!(h::Array{Float64,3}, L::Int, K::Float64, sigma::Float64;
                     therm_sweeps::Int=10000, use_sw::Bool=true, step_size::Float64=1.0)
    for i in 1:therm_sweeps
        metropolis_sweep!(h, L, K, sigma, step_size)
        if use_sw && i % 5 == 0
            swendsen_wang_sweep!(h, L, K, sigma)
        end
    end
end

end # module
