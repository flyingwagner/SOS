#!/bin/bash
#
# 并行运行所有参数组合，每个组合占用一个核
#

JULIA=~/.juliaup/bin/julia

L_LIST=(256 512)
SIGMA_LIST=(0.5 1.0 2.0 4.0)

for L in "${L_LIST[@]}"; do
    for sigma in "${SIGMA_LIST[@]}"; do
        echo "启动: L=$L, sigma=$sigma"
        $JULIA run_measure_extra.jl $L $sigma &
    done
done

echo "共启动 $((${#L_LIST[@]} * ${#SIGMA_LIST[@]})) 个并行任务，等待全部完成..."
wait
echo "全部完成"
