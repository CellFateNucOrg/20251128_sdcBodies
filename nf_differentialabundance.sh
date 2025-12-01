#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

genomeVer=WS295
genomeDir=/mnt/external.data/MeisterLab/publicData/genomes/${genomeVer}
gtfFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.canonical_geneset.protein_coding.gtf

WORK_DIR=/mnt/external.data/MeisterLab/jsemple/20251118_sdcBodies
RESULTS_DIR=$WORK_DIR/res18_no1298IPB2_LowInput_minAbun10_minSamples16_pc_noOsc_noShrink
mkdir -p $RESULTS_DIR
CONFIG_FILE=/mnt/external.data/MeisterLab/nf-core/unibe_izb.config

nextflow run nf-core/differentialabundance -profile rnaseq,singularity \
	--input ${WORK_DIR}/fileList_extended_compareLI.csv --outdir $RESULTS_DIR -c $CONFIG_FILE -r 1.5.0 \
	--gtf $gtfFile   \
	--matrix ${WORK_DIR}/star_salmon/salmon.merged.gene_counts.bulk_plus_lowinput.protein_coding_noOsc.tsv \
	--transcript_length_matrix ${WORK_DIR}/star_salmon/salmon.merged.gene_lengths.bulk_plus_lowinput.tsv \
	--contrasts ${WORK_DIR}/contrasts_compareLI.csv \
	--deseq2_shrink_lfc false  --filtering_min_abundance 10 --filtering_min_samples 16

#keep a record of this file in the results directory
cp nf_differentialabundance.sh $RESULTS_DIR/

