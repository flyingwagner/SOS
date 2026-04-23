#!/bin/bash
#SBATCH --job-name=SOS_extra
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --chdir=/home/shhu/sihan/SOS
#SBATCH --output=/home/shhu/sihan/SOS/logs/SOS_extra_%j.out
#SBATCH --error=/home/shhu/sihan/SOS/logs/SOS_extra_%j.err

mkdir -p /home/shhu/sihan/SOS/logs /home/shhu/sihan/SOS/data_extra

echo "任务开始: $(date)"
echo "节点: $(hostname)"
echo "分配CPU核数: $SLURM_CPUS_PER_TASK"

bash /home/shhu/sihan/SOS/run_extra_all.sh

echo "全部任务结束: $(date)"
