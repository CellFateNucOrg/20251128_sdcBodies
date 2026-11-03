#!/bin/bash
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1



WORK_DIR=/mnt/meister.data/sinem/sequencing/20261028_se_ill_rna_l3_wPM1556_autosomal_bmm1
RESULTS_DIR=$WORK_DIR/diffabund

pixi run summarise --params ../params.yml ../params_raptor.yml --workdir $WORK_DIR -o $RESULTS_DIR/da_summary.pdf
