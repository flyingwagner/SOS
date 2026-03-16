module BlockStats

using Statistics

export block_stats

"""
    block_stats(series, block_size; tol=0.1, min_blocks=64) -> (mean, err, rho1, final_block_size)

对时间序列做 blocking 分析，带自适应 coarsening。

若相邻 block 自相关 ρ(1) > tol，则将相邻两个 block 合并（block 大小翻倍），
重复直到 ρ(1) ≤ tol 或 block 数降至 min_blocks 以下。

返回：
  mean             : block 均值的均值
  err              : 标准误 std/√n_blocks
  rho1             : 最终 block 序列的相邻自相关系数
  final_block_size : 最终使用的 block 大小
"""
function block_stats(series::Vector{Float64}, block_size::Int;
                     tol::Float64=0.1, min_blocks::Int=64)
    n = div(length(series), block_size)
    # 初始 block 均值
    blocks = [mean(series[(b-1)*block_size+1 : b*block_size]) for b in 1:n]
    cur_block_size = block_size

    while true
        n_b = length(blocks)
        m   = mean(blocks)

        # 计算方差和 lag-1 协方差（滚动，仿 statistics.c）
        var_sum = 0.0
        cor_sum = 0.0
        devp    = 0.0
        for i in 1:n_b
            devn     = blocks[i] - m
            var_sum += devn * devn
            cor_sum += devn * devp
            devp     = devn
        end
        rho1 = (var_sum > 0) ? cor_sum / var_sum : 0.0

        if rho1 <= tol || n_b <= min_blocks
            err = sqrt(var_sum / (n_b - 1)) / sqrt(n_b)
            return m, err, rho1, cur_block_size
        end

        # Coarsen：合并相邻两个 block
        n_new = div(n_b, 2)
        blocks = [(blocks[2i-1] + blocks[2i]) / 2.0 for i in 1:n_new]
        cur_block_size *= 2
    end
end

end # module
