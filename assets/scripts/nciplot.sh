#!/usr/bin/env bash
sbatch <<EOT
#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH -t 0-24:00:00
#SBATCH -J NCIPLOT_$1
#SBATCH -o log.%j
# SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo '-------------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo Input file is $1, job name is $1
echo '-------------------------------'
echo -e '\n\n'

export PROCS=${SLURM_CPUS_ON_NODE}
export OMP_NUM_THREADS=32

/home/nt383/project/installers/src_nciplot_4.2.1_alpha/nciplot $1.nci

EOT
