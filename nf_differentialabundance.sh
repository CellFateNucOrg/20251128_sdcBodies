#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

genomeVer=WS298
genomeDir=/mnt/meister.data/publicData/genomes/${genomeVer}
gtfFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.canonical_geneset.gtf

WORK_DIR=/mnt/meister.data/sinem/sequencing/20261028_se_ill_rna_l3_wPM1556_autosomal_bmm1
RESULTS_DIR=$WORK_DIR/diffabund
mkdir -p $RESULTS_DIR
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb_apptainer.config

## differntial abundance runs with all samples, +-lfc, +-protein_coding +-oscillating
#nextflow run nf-core/differentialabundance -profile rnaseq,apptainer \
# 	--input ${WORK_DIR}/samplesheet.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 2.0.0 \
# 	--gtf $gtfFile   \
# 	--contrasts ${WORK_DIR}/contrasts.csv \
# 	--paramsheet params.yml \
# 	--paramset_name 01_noShrink_minAbund10_minSamples3,02_lfcShrink_minAbund10_minSamples3,03_noShrink_minAbund10_minSamples3_pc,04_lfcShrink_minAbund10_minSamples3_pc,05_noShrink_minAbund10_minSamples3_pc_noOsc,06_lfcShrink_minAbund10_minSamples3_pc_noOsc


## differntial abundance runs with all samples, +-lfc, +-protein_coding +-oscillating using raptor age estimations
nextflow run nf-core/differentialabundance -profile rnaseq,apptainer \
 	--input ${WORK_DIR}/samplesheet_raptor.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 2.0.0 \
 	--gtf $gtfFile   \
 	--contrasts ${WORK_DIR}/contrasts_raptor.yaml \
 	--paramsheet params_raptor.yml \
 	--paramset_name 01_noShrink_minAbund10_minSamples3_raptor,02_lfcShrink_minAbund10_minSamples3_raptor,03_noShrink_minAbund10_minSamples3_pc_raptor,04_lfcShrink_minAbund10_minSamples3_pc_raptor 


## differntial abundance runs with all samples, +-lfc, +-protein_coding +-oscillating using raptor age estimations after removing outlier samples
removeOutliers=false
if $removeOutliers; then
  nextflow run nf-core/differentialabundance -profile rnaseq,apptainer \
         --input ${WORK_DIR}/samplesheet_noOut.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 2.0.0 \
         --gtf $gtfFile   \
         --contrasts ${WORK_DIR}/contrasts.csv \
         --paramsheet params_noOut.yml \
 	--paramset_name 01_noShrink_minAbund10_minSamples3_noOut,02_lfcShrink_minAbund10_minSamples3_noOut,03_noShrink_minAbund10_minSamples3_pc_noOut,04_lfcShrink_minAbund10_minSamples3_pc_noOut,05_noShrink_minAbund10_minSamples3_pc_noOsc_noOut,06_lfcShrink_minAbund10_minSamples3_pc_noOsc_noOut

 

## differntial abundance runs with all samples, +-lfc, +-protein_coding +-oscillating using raptor age estimations after removing outlier samples
  nextflow run nf-core/differentialabundance -profile rnaseq,apptainer \
        --input ${WORK_DIR}/samplesheet_raptor_noOut.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 2.0.0 \
        --gtf $gtfFile   \
        --contrasts ${WORK_DIR}/contrasts_raptor.yaml \
        --paramsheet params_raptor.yml \
        --paramset_name 01_noShrink_minAbund10_minSamples3_raptor
fi


#keep a record of this file in the results directory
cp nf_differentialabundance.sh $RESULTS_DIR/
