module GeneralizedSOS

using Random
using Statistics
using PyCall
using PyPlot

# ==========================================
# 模拟参数
# ==========================================
const THERM_SWEEPS   = 10000   # 热化步数
const MEASURE_SWEEPS = 100000  # 测量步数

# ==========================================
# 模块 1：系统定义与核心更新算法
# ==========================================

# 简易高效的并查集 (Union-Find)，用于 Swendsen-Wang 集团算法
mutable struct UnionFind
    parent::Vector{Int}
    UnionFind(n::Int) = new(collect(1:n))
end

function find!(uf::UnionFind, i::Int)
    if uf.parent[i] == i
        return i
    end
    uf.parent[i] = find!(uf, uf.parent[i]) # 路径压缩
    return uf.parent[i]
end

function union!(uf::UnionFind, i::Int, j::Int)
    root_i = find!(uf, i)
    root_j = find!(uf, j)
    if root_i != root_j
        uf.parent[root_i] = root_j
    end
end

# 零模消除约束：强制整体质心高度为0
function zero_mode_shift!(h::Matrix{Float64})
    h .-= mean(h) #
end

# 获取周期性边界条件下的邻居索引
@inline function get_neighbors(x::Int, y::Int, L::Int)
    up    = y == L ? 1 : y + 1
    down  = y == 1 ? L : y - 1
    right = x == L ? 1 : x + 1
    left  = x == 1 ? L : x - 1
    return (x, up), (x, down), (right, y), (left, y)
end

# 单步 Metropolis 局部更新
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
            delta_E += K * (abs(new_h - neighbor_h)^sigma - abs(old_h - neighbor_h)^sigma) #
        end
        
        if delta_E <= 0.0 || rand() < exp(-delta_E)
            h[x, y] = new_h
            accepted += 1
        end
    end
    zero_mode_shift!(h) # 演化阶段定期平移防止数值溢出
    return accepted / N
end

# 单步 Swendsen-Wang 嵌入式集团更新 (Embedded Cluster)
# 随机选取反射平面 h_ref，将连续高度场投影为 Ising 自旋后做 FK 分解
function swendsen_wang_sweep!(h::Matrix{Float64}, L::Int, K::Float64, sigma::Float64)
    N = L * L
    uf = UnionFind(N)

    # 随机选取反射平面：取一个随机格点的高度
    h_ref = h[rand(1:L), rand(1:L)]

    # 将二维坐标映射到一维索引
    idx(x, y) = (y - 1) * L + x

    # 投影自旋 s_i = sign(h_i - h_ref)
    # 反射操作: h_i -> 2*h_ref - h_i
    # 只有同侧（同号投影自旋）的邻居才可能放键
    for y in 1:L
        for x in 1:L
            current_idx = idx(x, y)
            hi = h[x, y]
            si = hi > h_ref ? 1 : -1

            # 右侧邻居
            nx = x == L ? 1 : x + 1
            hj = h[nx, y]
            sj = hj > h_ref ? 1 : -1
            if si == sj  # 同侧才放键
                # 反射前能量: K|hi - hj|^sigma
                # 反射 i 后能量: K|(2*h_ref - hi) - hj|^sigma
                E_before = K * abs(hi - hj)^sigma
                E_after  = K * abs(2.0 * h_ref - hi - hj)^sigma
                delta_E = E_after - E_before
                if delta_E > 0.0
                    p_bond = 1.0 - exp(-delta_E)
                    if rand() < p_bond
                        union!(uf, current_idx, idx(nx, y))
                    end
                end
            end

            # 上方邻居
            ny = y == L ? 1 : y + 1
            hj = h[x, ny]
            sj = hj > h_ref ? 1 : -1
            if si == sj
                E_before = K * abs(hi - hj)^sigma
                E_after  = K * abs(2.0 * h_ref - hi - hj)^sigma
                delta_E = E_after - E_before
                if delta_E > 0.0
                    p_bond = 1.0 - exp(-delta_E)
                    if rand() < p_bond
                        union!(uf, current_idx, idx(x, ny))
                    end
                end
            end
        end
    end

    # 为每个集团以 50% 概率做反射 h_i -> 2*h_ref - h_i
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
# 模块 2：物理量测量模块
# ==========================================

# 测量能量密度
function measure_energy(h::Matrix{Float64}, L::Int, K::Float64, sigma::Float64)
    E = 0.0
    for y in 1:L
        for x in 1:L
            nx = x == L ? 1 : x + 1
            ny = y == L ? 1 : y + 1
            E += K * abs(h[x, y] - h[nx, y])^sigma # x方向键
            E += K * abs(h[x, y] - h[x, ny])^sigma # y方向键
        end
    end
    return E / (L * L)
end

# 测量表面粗糙度 W^2
function measure_roughness(h::Matrix{Float64})
    return mean(h.^2) # 均值已通过 zero_mode_shift 保证为0
end

# 测量沿 x 方向的高度-高度关联函数 G(r)
function measure_correlation(h::Matrix{Float64}, L::Int)
    max_r = div(L, 2)
    G = zeros(Float64, max_r)
    for r in 1:max_r
        sum_diff = 0.0
        for y in 1:L
            for x in 1:L
                nx = mod1(x + r, L)
                sum_diff += (h[nx, y] - h[x, y])^2 #
            end
        end
        G[r] = sum_diff / (L * L)
    end
    return G
end

# 测量动量空间的结构因子 chi_k
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
    
    # 按照定义 |h_k|^2 返回
    return (sum_cos^2 + sum_sin^2) / (L * L)
end


# ==========================================
# 模块 3：主模拟流程与控制器
# ==========================================

function run_simulation(; L::Int=32, sigma::Float64=2.0, K::Float64=1.0,
                          therm_sweeps::Int=THERM_SWEEPS, measure_sweeps::Int=MEASURE_SWEEPS,
                          use_sw::Bool=true)
    h = zeros(Float64, L, L) # 初始化平坦表面
    step_size = 1.0 # Metropolis 步长，可自适应调节
    
    # 物理量存储器
    energies = Float64[]
    roughnesses = Float64[]
    G_avg = zeros(Float64, div(L, 2))
    chi_kmin_avg = 0.0
    
    println("--- 开始模拟 L=$L, sigma=$sigma, 混合算法=$use_sw ---")
    
    # 1. 热化阶段 (Thermalization)
    for i in 1:therm_sweeps
        metropolis_sweep!(h, L, K, sigma, step_size)
        if use_sw && i % 5 == 0
            swendsen_wang_sweep!(h, L, K, sigma)
        end
    end
    
    # 2. 测量阶段 (Measurement)
    for i in 1:measure_sweeps
        metropolis_sweep!(h, L, K, sigma, step_size)
        if use_sw && i % 5 == 0
            swendsen_wang_sweep!(h, L, K, sigma)
        end
        
        # 测量采样
        push!(energies, measure_energy(h, L, K, sigma))
        push!(roughnesses, measure_roughness(h))
        G_avg .+= measure_correlation(h, L)
        chi_kmin_avg += measure_structure_factor_kmin(h, L)
    end
    
    # 平均化
    G_avg ./= measure_sweeps
    chi_kmin_avg /= measure_sweeps
    
    # 计算比热 Cv = (<H^2> - <H>^2) / (N * T^2)  (这里 T 吸收到 K 中，设为1)
    Cv = L * L * var(energies) #
    
    return mean(energies), Cv, mean(roughnesses), G_avg, chi_kmin_avg
end


# ==========================================
# 模块 4：科学作图与有限尺寸标度分析
# ==========================================

function run_and_plot_results()
    L_list = [8, 16, 32] 
    sigma_list = [1.0, 2.0, 3.0] 
    
    # 存储粗糙度数据
    W2_results = Dict{Float64, Vector{Float64}}()
    
    for sigma in sigma_list
        W2_results[sigma] = Float64[]
        for L in L_list
            _, _, w2, _, _ = run_simulation(L=L, sigma=sigma)
            push!(W2_results[sigma], w2)
        end
    end
    
    # -----------------------------------------------------
    # 引入 Python 的 scipy 进行标度拟合 (Scaling Fit)
    # -----------------------------------------------------
    scipy_stats = pyimport("scipy.stats")
    println("\n--- 有限尺寸标度 (FSS) 拟合结果 ---")
    
    # 画图设置初始化 (使用 matplotlib.pyplot API)
    fig = figure("SOS Scaling Results", figsize=(10, 4))
    
    # 子图 1: 粗糙度 W^2 随尺寸 L 的变化
    subplot(1, 2, 1)
    title("Surface Roughness Scaling")
    xlabel("System Size L (log scale)")
    ylabel("\$W^2(L)\$")
    xscale("log", base=10) # 设置对数坐标
    
    for sigma in sigma_list
        # PyPlot 画线
        plot(L_list, W2_results[sigma], marker="o", label="\$\\sigma = $sigma\$", linewidth=2)
        log_L = log.(L_list)
        # 使用 scipy.stats.linregress 对 W^2 和 ln(L) 进行线性拟合
        slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(log_L, W2_results[sigma])
        println("\$\\sigma\$ = $sigma 标度拟合: W² ~ $(round(slope, digits=3)) * ln(L) + $(round(intercept, digits=3))")
    end
    legend(loc="upper left")
    
    # 子图 2: 选取 L=32 观察 G(r)
    L_target = 32
    subplot(1, 2, 2)
    title("Height-Height Correlation G(r) for L=$L_target")
    xlabel("Distance r")
    ylabel("G(r)")
    
    for sigma in sigma_list
        _, _, _, G, _ = run_simulation(L=L_target, sigma=sigma)
        plot(1:div(L_target, 2), G, marker="s", label="\$\\sigma = $sigma\$", linewidth=2)
    end
    legend(loc="upper left")
    
    # 调整布局并保存
    tight_layout()
    savefig("SOS_Scaling_Results.png", dpi=300)
    savefig("SOS_Scaling_Results.pdf")
    println("\n图表已使用 PyPlot 保存为 SOS_Scaling_Results.png 和 .pdf")
end

end # module 结束

if abspath(PROGRAM_FILE) == @__FILE__
    using .GeneralizedSOS
    
    # 2. 加上模块名前缀来调用函数
    GeneralizedSOS.run_and_plot_results()
end