#!/bin/bash
#
# 并行运行所有 3D (L, σ) 组合，每个组合占用一个核
#

JULIA=~/.juliaup/bin/julia
WORKDIR=/home/shhu/sihan/SOS
DATADIR=$WORKDIR/data3d

# 单进程内不让 BLAS / OpenMP 抢核，避免 16 个 Julia 进程互相挤
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_NUM_THREADS=1

L_LIST=(4 8 16 32)
SIGMA_LIST=(0.5 1.0 2.0 4.0)

for L in "${L_LIST[@]}"; do
    for sigma in "${SIGMA_LIST[@]}"; do
        echo "启动: L=$L, sigma=$sigma"
        $JULIA $WORKDIR/run3d.jl $L $sigma --datadir $DATADIR \
            > $WORKDIR/logs/run3d_L${L}_s${sigma}.log 2>&1 &
    done
done

echo "共启动 $((${#L_LIST[@]} * ${#SIGMA_LIST[@]})) 个并行任务，等待全部完成..."
wait
echo "全部完成"
