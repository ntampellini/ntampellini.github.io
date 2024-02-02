#!/bin/bash
sbatch <<EOT
#!/usr/bin/env bash

# SBATCH --ntasks=$2
# SBATCH --cpus-per-task=1

#SBATCH --nodes=1
#SBATCH --tasks-per-node=$2
#SBATCH --mem=$((5000*$2))
#SBATCH -t 0-24:00:00
#SBATCH -J Orca_${PWD##*/}_$1
#SBATCH -o log.%j
# SBATCH --mail-type=begin
#SBATCH --mail-type=end


echo '-------------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo Job name is $1
echo '-------------------------------'
echo -e '\n\n'

module purge
module load xtb/6.5.1-foss-2020b ORCA/5.0.4-gompi-2020b
module list

export RSH_COMMAND="/usr/bin/ssh -x"
export LD_LIBRARY_PATH=/vast/palmer/apps/avx2/software/GCCcore/10.2.0/lib64/:/vast/palmer/apps/avx2/software/ORCA/5.0.4-gompi-2020b/lib/:/vast/palmer/apps/avx2/software/ORCA/5.0.4-gompi-2020b/bin/
/vast/palmer/apps/avx2/software/ORCA/5.0.4-gompi-2020b/bin/orca $1.inp > $1.out

EOT