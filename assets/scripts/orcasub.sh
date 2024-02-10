#!/bin/bash

# set options
gb_per_core=8

############################################### WRAPPER

# Check the number of arguments
if [ "$#" -lt 1 ]; then
    echo "  Not enough arguments! At least one is required."
    echo "  Usage: $0 input[.inp/.xyz] [n_cores]"
    echo "  If ncores is not specified, it is tentatively extracted from input.inp"
    exit 1
fi

# remove extension from base name
basename="${1%.*}"

# make sure input file
if ! [ -e "$basename.inp" ]; then
    echo "File $basename.inp not found. Exiting."
    exit 1
fi

# if only one argument was provided, try to extract the number
# of processors from the input file
if [ "$#" -lt 2 ]; then
    nprocs=$(awk '/nprocs/ {print $2}' "$basename.inp")
else 
    nprocs=$2
fi

# check that variable nprocs has been
# defined in a way or another
if ! [[ -v nprocs ]]; then
    echo "Failed extracting number of cores from input file. Please specify it manually."
    exit 1
fi

memory=$(($gb_per_core*$nprocs))G

echo "Running $basename on $nprocs CPUs/$memory ($gb_per_core GB/core)"

############################################### WRAPPED SCRIPT

sbatch <<EOT
#!/usr/bin/env bash

# SBATCH --ntasks=$nprocs
# SBATCH --cpus-per-task=1

#SBATCH --nodes=1
#SBATCH --tasks-per-node=$nprocs
#SBATCH --mem=$memory
#SBATCH -t 0-24:00:00
#SBATCH -J Orca_${PWD##*/}_$basename
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

orcadir=/vast/palmer/apps/avx2/software/ORCA/5.0.4-gompi-2020b/bin/

export RSH_COMMAND="/usr/bin/ssh -x"
export LD_LIBRARY_PATH=/vast/palmer/apps/avx2/software/GCCcore/10.2.0/lib64/:/vast/palmer/apps/avx2/software/ORCA/5.0.4-gompi-2020b/lib/:$orcadir

$orcadir/orca "${basename}.inp" > "${basename}.out"

EOT