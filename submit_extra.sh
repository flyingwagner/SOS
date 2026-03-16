#!/bin/bash
#
# 批量提交 run_measure_extra.jl 任务
# 参数组合: L = 256, 512; sigma = 0.5, 1.0, 2.0, 4.0
#

mkdir -p logs data_extra

L_LIST=(256 512)
SIGMA_LIST=(0.5 1.0 2.0 4.0)

for L in "${L_LIST[@]}"; do
    for sigma in "${SIGMA_LIST[@]}"; do
        job_name="SOS_L${L}_s${sigma}"
        echo "提交任务: $job_name"

        sbatch --job-name="$job_name" \
               --nodes=1 \
               --ntasks-per-node=1 \
               --cpus-per-task=1 \
               --time=24:00:00 \
               --output="logs/${job_name}_%j.out" \
               --error="logs/${job_name}_%j.err" \
               --wrap="julia run_measure_extra.jl $L $sigma"
    done
done

echo ""
echo "共提交 $((${#L_LIST[@]} * ${#SIGMA_LIST[@]})) 个任务"
echo "查看状态: squeue -u \$(whoami)"
