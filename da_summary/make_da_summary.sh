#!/bin/bash
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1



WORK_DIR=/mnt/meister.data/kalyan/sequencing/20261028_kg_ill_rna_l3_wPM1486_autosomal_smc1
RESULTS_DIR=$WORK_DIR/diffabund


PARAMFILES=( $WORK_DIR/params.yml $WORK_DIR/params_raptor.yml )
removeOutliers=true
if $removeOutliers; then
 PARAMFILES=( ../params.yml ../params_raptor.yml ../params_noOut.yml ../params_raptor_noOut.yml)
fi

pixi run summarise --params ${PARAMFILES[@]}  --workdir $WORK_DIR -o $RESULTS_DIR/da_summary.pdf
