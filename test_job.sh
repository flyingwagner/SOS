#!/bin/bash
#SBATCH --job-name=SOS_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=/home/shhu/sihan/SOS/logs/test_%j.out
#SBATCH --error=/home/shhu/sihan/SOS/logs/test_%j.err
#SBATCH --chdir=/home/shhu/sihan/SOS

alias julia="~/.juliaup/bin/julia"
alias mpiexecjl="~/.julia/bin/mpiexecjl"

mkdir -p logs data

echo "测试任务开始: $(date)"
echo "节点: $(hostname)"

julia run.jl 16 2.0 --therm 1000 --measure 10000

echo "测试任务结束: $(date)"
