#!/usr/bin/env bash
sbatch <<EOT
#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=2G
#SBATCH -t 0-24:00:00
#SBATCH -J xtbopt_$2
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

export OMP_NUM_THREADS=${SLURM_CPUS_ON_NODE}
export MKL_NUM_THREADS=${SLURM_CPUS_ON_NODE}

echo 'Command used is python /gpfs/gibbs/project/miller/nt383/scripts/xtbopt.py $1 constr,$3 newfile > xtbopt_shell.out'

python /gpfs/gibbs/project/miller/nt383/scripts/xtbopt.py $1 constr,$3 newfile > xtbopt_shell.out

EOT
