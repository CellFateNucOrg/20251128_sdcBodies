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

WORK_DIR=/mnt/meister.data/sinem/sequencing/20251118_Sinem_sdcBodies
RESULTS_DIR=$WORK_DIR/diffabund
mkdir -p $RESULTS_DIR
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb_apptainer.config



PARAMSETS=$(sed -n 's/^[[:space:]-]*paramset_name:[[:space:]]*//p' params.yml | tr -d "\"'" | paste -sd, -)

# differntial abundance runs with all samples, +-lfc, +-protein_coding +-oscillating
nextflow run nf-core/differentialabundance -profile rnaseq,apptainer \
 	--input ${WORK_DIR}/samplesheet.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 2.0.0 \
 	--gtf $gtfFile   \
 	--contrasts ${WORK_DIR}/contrasts.csv \
 	--paramsheet params.yml \
 	--paramset_name "$PARAMSETS" -resume

PARAMSETS=$(sed -n 's/^[[:space:]-]*paramset_name:[[:space:]]*//p' params_raptor.yml | tr -d "\"'" | paste -sd, -)

## differntial abundance runs with all samples, +-lfc, +-protein_coding +-oscillating using raptor age estimations
nextflow run nf-core/differentialabundance -profile rnaseq,apptainer \
 	--input ${WORK_DIR}/samplesheet_raptor.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 2.0.0 \
 	--gtf $gtfFile   \
 	--contrasts ${WORK_DIR}/contrasts_raptor.yml \
 	--paramsheet params_raptor.yml \
 	--paramset_name "$PARAMSETS" -resume


## differntial abundance runs with all samples, +-lfc, +-protein_coding +-oscillating using raptor age estimations after removing outlier samples
removeOutliers=false
if $removeOutliers; then
  PARAMSETS=$(sed -n 's/^[[:space:]-]*paramset_name:[[:space:]]*//p' params_noOut.yml | tr -d "\"'" | paste -sd, -)
  nextflow run nf-core/differentialabundance -profile rnaseq,apptainer \
         --input ${WORK_DIR}/samplesheet_noOut.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 2.0.0 \
         --gtf $gtfFile   \
         --contrasts ${WORK_DIR}/contrasts.csv \
         --paramsheet params_noOut.yml \
 	--paramset_name "$PARAMSETS"


## differntial abundance runs with all samples, +-lfc, +-protein_coding +-oscillating using raptor age estimations after removing outlier samples
  PARAMSETS=$(sed -n 's/^[[:space:]-]*paramset_name:[[:space:]]*//p' params_raptor_noOut.yml | tr -d "\"'" | paste -sd, -)
  nextflow run nf-core/differentialabundance -profile rnaseq,apptainer \
        --input ${WORK_DIR}/samplesheet_raptor_noOut.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 2.0.0 \
        --gtf $gtfFile   \
        --contrasts ${WORK_DIR}/contrasts_raptor.yml \
        --paramsheet params_raptor_noOut.yml \
        --paramset_name "$PARAMSETS"
fi


#keep a record of this file in the results directory
cp nf_differentialabundance.sh $RESULTS_DIR/
cp params.yml params_raptor.yml ${WORK_DIR}/contrasts.csv ${WORK_DIR}/contrasts_raptor.yml $RESULTS_DIR/

