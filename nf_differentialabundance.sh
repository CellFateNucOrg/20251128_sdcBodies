#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

WORK_DIR=/mnt/external.data/MeisterLab/jsemple/20251118_sdcBodies


if [ -f "${WORK_DIR}/variables.R" ]
then
  echo "Using variables.R"
  FILE_VAR=($(Rscript ${WORK_DIR}/queryVariables.R))
  genomeVer=${FILE_VAR[1]#\"}
  gtfFile=${FILE_VAR[2]}
  RUN_NAME=${FILE_VAR[3]}
  SAMPLE_SHEET=${FILE_VAR[4]}
  COUNTS=${FILE_VAR[5]}
  LENGTHS=${FILE_VAR[6]}
  CONTRASTS=${FILE_VAR[7]}
  SHRINK=${FILE_VAR[8]}
  MIN_ABUND=${FILE_VAR[9]}
  MIN_SAMPLES=${FILE_VAR[10]%\"}
else
  echo "No variables file found, using defaults"
  genomeVer=WS295
  GENOME_DIR=/mnt/external.data/MeisterLab/publicData/genomes/${genomeVer}
  gtfFile=${GENOME_DIR}/c_elegans.PRJNA13758.${genomeVer}.canonical_geneset.protein_coding.noMito.noOsc.gtf
  SAMPLE_SHEET=${WORK_DIR}/fileList_extended_compareLI.csv
  COUNTS=${WORK_DIR}/star_salmon/salmon.merged.gene_counts.bulk_plus_lowinput.protein_coding_noOsc.tsv
  LENGTHS=${WORK_DIR}/star_salmon/salmon.merged.gene_lengths.bulk_plus_lowInput.tsv
  CONTRASTS=${WORK_DIR}/contrasts_compareLI.csv
  SHRINK=false
  MIN_ABUND=10
  MIN_SAMPLES=16
  RUN_NAME=res01_no1298IPB2_LowInput_minAbun${MIN_ABUND}_minSamples${MIN_SAMPLES}_pc_noOsc_noShrink
fi

echo "running pipeline with following variables:"
echo $gtfFile
echo $RUN_NAME
echo $SAMPLE_SHEET
echo $COUNTS
echo $LENGTHS
echo $CONTRASTS
echo $SHRINK
echo $MIN_ABUND
echo $MIN_SAMPLES

RESULTS_DIR=$WORK_DIR/$RUN_NAME
mkdir -p $RESULTS_DIR
CONFIG_FILE=/mnt/external.data/MeisterLab/nf-core/unibe_izb.config

Rscript --vanilla ./filterCountsMatrix.R

nextflow run nf-core/differentialabundance -profile rnaseq,singularity \
	--input $SAMPLE_SHEET --outdir $RESULTS_DIR -c $CONFIG_FILE -r 1.5.0 \
	--gtf $gtfFile   \
	--matrix  $COUNTS \
	--transcript_length_matrix $LENGTHS \
	--contrasts $CONTRASTS  \
	--deseq2_shrink_lfc $SHRINK  --filtering_min_abundance $MIN_ABUND --filtering_min_samples $MIN_SAMPLES 

#keep a record of this file in the results directory
cp nf_differentialabundance.sh $RESULTS_DIR/
if [ -f "${WORK_DIR}/variables.R" ]; then
  cp ${WORK_DIR}/variables.R $RESULTS_DIR/
fi

