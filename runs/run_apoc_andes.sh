#!/bin/bash
#SBATCH -A bip198
#SBATCH -p batch
#SBATCH -N 60
#SBATCH -J apoc
#SBATCH --mem=0
##SBATCH --no-requeue  # prevent automatic restart due to node failure
#SBATCH -t 32:00:00
#SBATCH -o logs/%x.%j.stdout
#SBATCH -e logs/%x.%j.stderr

module load intel
#module load pgi

#target_lst=$SLURM_SUBMIT_DIR/job/R_rubrum_lt1480aa_hyp.lst
target_lst=$SLURM_SUBMIT_DIR/job/pseudodesulfovibrio_run3.lst

NUM_RUNS=$SLURM_JOB_NUM_NODES
NUM_CORES_PER_RUN=$SLURM_CPUS_ON_NODE  # 32 cores per node at Andes

dat_root=/gpfs/alpine/world-shared/bif135/species/pseudodesulfovibrio/af_mod_genome
out_root=/gpfs/alpine/world-shared/bif135/species/temp/pseudodesulfovibrio/apoc_rbd

#export OMP_NUM_THREADS=32

cd $SLURM_SUBMIT_DIR || { echo "Error: could not reach the working directory" ; exit 1; }

# used to create sublists of query structures
../scripts/splitset.pl -s $target_lst $NUM_RUNS

start_time="$(date -u +%s)"
date

for (( i=1; i<=$NUM_RUNS; i++ ))
do
  echo "Info: run $i at $HOSTNAME"
  srun -N1 -n1 -c $NUM_CORES_PER_RUN ./apoc_pdb70_andes.pl -c $NUM_CORES_PER_RUN \
    -dat_root $dat_root \
    -out_root $out_root \
    ${target_lst}.${i} &
done
wait

date
end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Info: all runs are done, run time $elapsed seconds"

