module SOSCore

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
# 基础工具
# ==========================================
function zero_mode_shift!(h::Matrix{Float64})
    h .-= mean(h)
end

@inline function get_neighbors(x::Int, y::Int, L::Int)
    up    = y == L ? 1 : y + 1
    down  = y == 1 ? L : y - 1
    right = x == L ? 1 : x + 1
    left  = x == 1 ? L : x - 1
    return (x, up), (x, down), (right, y), (left, y)
end

# ==========================================
# 更新算法
# ==========================================
function metropolis_sweep!(h::Matrix{Float64}, L::Int, K::Float64, sigma::Float64, step_size::Float64)
    accepted = 0
    N = L * L
    for _ in 1:N
        x, y = rand(1:L), rand(1:L)
        old_h = h[x, y]
        new_h = old_h + step_size * (2.0 * rand() - 1.0)

        delta_E = 0.0
        for (nx, ny) in get_neighbors(x, y, L)
            neighbor_h = h[nx, ny]
            delta_E += K * (abs(new_h - neighbor_h)^sigma - abs(old_h - neighbor_h)^sigma)
        end

        if delta_E <= 0.0 || rand() < exp(-delta_E)
            h[x, y] = new_h
            accepted += 1
        end
    end
    zero_mode_shift!(h)
    return accepted / N
end

function swendsen_wang_sweep!(h::Matrix{Float64}, L::Int, K::Float64, sigma::Float64)
    N = L * L
    uf = UnionFind(N)
    h_ref = h[rand(1:L), rand(1:L)]
    idx(x, y) = (y - 1) * L + x

    for y in 1:L
        for x in 1:L
            current_idx = idx(x, y)
            hi = h[x, y]
            si = hi > h_ref ? 1 : -1

            nx = x == L ? 1 : x + 1
            hj = h[nx, y]
            sj = hj > h_ref ? 1 : -1
            if si == sj
                delta_E = K * abs(2.0 * h_ref - hi - hj)^sigma - K * abs(hi - hj)^sigma
                if delta_E > 0.0 && rand() < 1.0 - exp(-delta_E)
                    union!(uf, current_idx, idx(nx, y))
                end
            end

            ny = y == L ? 1 : y + 1
            hj = h[x, ny]
            sj = hj > h_ref ? 1 : -1
            if si == sj
                delta_E = K * abs(2.0 * h_ref - hi - hj)^sigma - K * abs(hi - hj)^sigma
                if delta_E > 0.0 && rand() < 1.0 - exp(-delta_E)
                    union!(uf, current_idx, idx(x, ny))
                end
            end
        end
    end

    cluster_flip = Dict{Int, Bool}()
    for i in 1:N
        root = find!(uf, i)
        if !haskey(cluster_flip, root)
            cluster_flip[root] = rand() < 0.5
        end
        if cluster_flip[root]
            cy = div(i - 1, L) + 1
            cx = mod(i - 1, L) + 1
            h[cx, cy] = 2.0 * h_ref - h[cx, cy]
        end
    end
    zero_mode_shift!(h)
end

# ==========================================
# 热化 + 演化的通用 sweep
# ==========================================
"""
    thermalize!(h, L, K, sigma; therm_sweeps, use_sw, step_size)

热化阶段：执行 therm_sweeps 步 Metropolis，每 5 步插入一次 SW。
"""
function thermalize!(h::Matrix{Float64}, L::Int, K::Float64, sigma::Float64;
                     therm_sweeps::Int=10000, use_sw::Bool=true, step_size::Float64=1.0)
    for i in 1:therm_sweeps
        metropolis_sweep!(h, L, K, sigma, step_size)
        if use_sw && i % 5 == 0
            swendsen_wang_sweep!(h, L, K, sigma)
        end
    end
end

end # module
