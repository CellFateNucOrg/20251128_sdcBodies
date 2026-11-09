#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

genomeVer=WS298
genomeDir=/mnt/meister.data/publicData/genomes/${genomeVer}
genomeFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.genomic.fa
gtfFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.canonical_geneset.gtf

WORK_DIR=/mnt/meister.data/sinem/sequencing/20251118_Sinem_sdcBodies/
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb_apptainer.config


nextflow run nf-core/rnaseq -profile apptainer -r 3.27.0 --input ${WORK_DIR}/samplesheet.csv --multiqc_title multiqc_rnaseq --outdir $WORK_DIR -c $CONFIG_FILE --fasta $genomeFile --gtf $gtfFile 

