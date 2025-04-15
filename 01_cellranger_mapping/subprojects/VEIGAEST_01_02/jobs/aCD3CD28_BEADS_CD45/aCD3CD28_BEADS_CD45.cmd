#!/usr/bin/env bash

#SBATCH --job-name="aCD3CD28_BEADS_CD45"
#SBATCH --chdir=.

#SBATCH --error=./logs/slurm_%x_%J.err
#SBATCH --output=./logs/slurm_%x_%J.out

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=marathon
#SBATCH --time=6-11:59:00
#SBATCH --partition=genD
#SBATCH --mem=2G

echo [`date "+%Y-%m-%d %T"`] started job on $HOSTNAME

export HDF5_USE_FILE_LOCKING="FALSE"

/scratch_isilon/groups/singlecell/shared/software/cellranger/cellranger-9.0.0/cellranger multi --id aCD3CD28_BEADS_CD45 \
    --csv /scratch_isilon/groups/singlecell/ljimenez/projects/VEIGAEST/01_cellranger_mapping/subprojects/VEIGAEST_01_02/jobs/aCD3CD28_BEADS_CD45/aCD3CD28_BEADS_CD45_config.csv \
    --localcores 2 \
    --jobmode /scratch_isilon/groups/singlecell/shared/software/cellranger/cellranger-9.0.0/external/martian/jobmanagers/slurm.template

echo [`date "+%Y-%m-%d %T"`] finished job
