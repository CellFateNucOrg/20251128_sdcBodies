#!/bin/bash
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1



WORK_DIR=/mnt/meister.data/sinem/sequencing/20251118_Sinem_sdcBodies
RESULTS_DIR=$WORK_DIR/diffabund


PARAMFILES=( ../params.yml ../params_raptor.yml )
removeOutliers=false
if $removeOutliers; then
 PARAMFILES=( ../params.yml ../params_raptor.yml ../params_noOut.yml ../params_raptor_noOut.yml)
fi

pixi run summarise --params ${PARAMFILES[@]}  --workdir $WORK_DIR -o $RESULTS_DIR/da_summary.pdf
