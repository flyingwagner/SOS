#!/bin/bash
#
# L=128, 全部 σ 并行（4 个进程，每个一核）
#

JULIA=~/.juliaup/bin/julia
WORKDIR=/home/shhu/sihan/SOS
DATADIR=$WORKDIR/data3d

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_NUM_THREADS=1

L=128
SIGMA_LIST=(0.5 1.0 2.0 4.0)

for sigma in "${SIGMA_LIST[@]}"; do
    echo "启动: L=$L, sigma=$sigma"
    $JULIA $WORKDIR/run3d.jl $L $sigma --datadir $DATADIR \
        > $WORKDIR/logs/run3d_L${L}_s${sigma}.log 2>&1 &
done

echo "共启动 ${#SIGMA_LIST[@]} 个并行任务，等待全部完成..."
wait
echo "全部完成"
