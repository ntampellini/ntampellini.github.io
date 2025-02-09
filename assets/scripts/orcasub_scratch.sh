#!/bin/bash

#   sbatch <<EOT
#   !/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --mem=128G
#SBATCH -t 0-24:00:00
#SBATCH -J Orca_scratchtest
#SBATCH -o log.%j
# SBATCH --mail-type=begin
#SBATCH --mail-type=end

# tell Slurm to send the SIGUSR1 signal two minutes before timeout
#SBATCH --signal=B:SIGUSR1@120


# echo '-------------------------------'
# cd ${SLURM_SUBMIT_DIR}
# echo ${SLURM_SUBMIT_DIR}
# echo Running on host $(hostname)
# echo Time is $(date)
# echo SLURM_NODES are $(echo ${SLURM_NODELIST})
# echo Job name is $1
# echo '-------------------------------'
# echo -e '\n\n'

function copy_files_to_submit_dir {
  cp $tdir/$job.out $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/$job.gbw $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.engrad $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/$job*.xyz $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.loc $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.qro $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.uno $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.unso $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.uco $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/$job.hess $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.cis $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.dat $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.mp2nat $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.nat $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.scfp_fod $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.scfp $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.scfr $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.nbo $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/FILE.47 $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.txt $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*spin* $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.densities* $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.nto $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.dat $SLURM_SUBMIT_DIR 2>/dev/null || :
  cp $tdir/*.allxyz $SLURM_SUBMIT_DIR 2>/dev/null || :
}

trap '"Timeout! Transfering files back to ${SLURM_SUBMIT_DIR} before Slurm kills the job." >> ${1}.out; copy_files_to_submit_dir; exit 2' SIGUSR1

module purge
module load xtb/6.6.0-foss-2020b ORCA/6.0.0-gompi-2022b

export RSH_COMMAND="/usr/bin/ssh -x"
export orcadir=/vast/palmer/apps/avx2/software/ORCA/6.0.0-gompi-2022b/bin
export LD_LIBRARY_PATH=/vast/palmer/apps/avx2/software/GCCcore/10.2.0/lib64/:/vast/palmer/apps/avx2/software/ORCA/6.0.0-gompi-2022b/lib/:$orcadir
export job=$1

# let ORCA know where the XTB executable is
export XTBEXE=/vast/palmer/apps/avx2/software/xtb/6.6.0-foss-2020b/bin/xtb

# Creating local scratch folder for the user on the computing node, if none exists.
export scratchlocation=/vast/palmer/scratch/miller

if [ ! -d $scratchlocation/$USER ]
then
  mkdir -p $scratchlocation/$USER
fi

tdir=$(mktemp -d $scratchlocation/$USER/ORCA_${PWD##*/}-${job}__$SLURM_JOB_ID-XXXX)

# Copy only the necessary stuff in submit directory to scratch directory. Add more here if needed.
cp  $SLURM_SUBMIT_DIR/$job.inp $tdir/
cp  $SLURM_SUBMIT_DIR/$job.gbw $tdir/ 2>/dev/null || :
cp  $SLURM_SUBMIT_DIR/*.cmp $tdir/ 2>/dev/null ||:
cp  $SLURM_SUBMIT_DIR/*.oldgbw $tdir/ 2>/dev/null || :
cp  $SLURM_SUBMIT_DIR/*.hess $tdir/ 2>/dev/null || :

# transfers all .xyz files and .allxyz files, except trajectories - might be an overkill in most cases but it is required for restarting NEBs
# cp  $SLURM_SUBMIT_DIR/*[^trj].xyz $tdir/ 2>/dev/null || :
find $SLURM_SUBMIT_DIR/ -type f -name "*.xyz" ! -name "*trj*" -exec cp {} $tdir \; >/dev/null || :


# Creating nodefile in scratch
echo $SLURM_NODELIST > $tdir/$job.nodes

# cd to scratch
cd $tdir

# Copy job and node info to beginning of outputfile
echo "Job execution start: $(date)" >> $job.out
echo "Shared library path: $LD_LIBRARY_PATH" >> $job.out
echo "Slurm Job ID is: ${SLURM_JOB_ID}" >> $job.out
echo "Slurm Job name is: ${SLURM_JOB_NAME}" >> $job.out
echo $SLURM_NODELIST >> $job.out

# Start ORCA job. ORCA is started using full pathname (necessary for parallel execution).
# Keeping .out file in the scratch directory so that the job can be checked with check.py
# $orcadir/orca $job.inp >> $SLURM_SUBMIT_DIR/$job.out
$orcadir/orca $job.inp >> $job.out &

# wait for the background command to finish or the USRSIG1 signal
wait

# ORCA has finished here. Now copy important stuff back (xyz files, GBW files etc.).
copy_files_to_submit_dir

exit 0
