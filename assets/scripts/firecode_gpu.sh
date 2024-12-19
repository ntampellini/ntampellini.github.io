#!/bin/bash
sbatch <<EOT
#!/usr/bin/env bash

#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --mem-per-cpu=4000
#SBATCH -t 0-24:00:00
#SBATCH -J 🔥_FIRECODE_GPU_$2
#SBATCH -o log.%j
# SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo '-------------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo Input file is $1, job name is $2
echo '-------------------------------'
echo -e '\n\n'

export OMP_NUM_THREADS=4,1
export OMP_MAX_ACTIVE_LEVELS=1
export MKL_NUM_THREADS=4

python -m firecode $1 -n $2

EOT
