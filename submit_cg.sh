#!/bin/bash
#SBATCH --job-name=SOS_cg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --chdir=/home/shhu/sihan/SOS
#SBATCH --output=/home/shhu/sihan/SOS/logs/SOS_cg_%j.out
#SBATCH --error=/home/shhu/sihan/SOS/logs/SOS_cg_%j.err

mkdir -p /home/shhu/sihan/SOS/logs /home/shhu/sihan/SOS/data_cg

echo "任务开始: $(date)"
echo "节点: $(hostname)"
echo "分配CPU核数: $SLURM_CPUS_PER_TASK"

JULIA=~/.juliaup/bin/julia
SIGMA_LIST=(0.5 1.0 2.0 4.0)

for sigma in "${SIGMA_LIST[@]}"; do
    echo "启动: L=1024, sigma=$sigma"
    $JULIA /home/shhu/sihan/SOS/run_coarse_grain.jl 1024 $sigma \
        --datadir /home/shhu/sihan/SOS/data_cg &
done

echo "共启动 ${#SIGMA_LIST[@]} 个并行任务，等待全部完成..."
wait
echo "全部任务结束: $(date)"
