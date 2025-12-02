# 20251128_sdcBodies

## RNAseq alignment

Sample sheet and contrast matrix created with ***./makeSampleSheet.R***

Initial alignment to genome/transcriptome carried out with ***./nf_rnaseq_submit.sh***

## Differential expression analysis

Results of DE analysis are each placed into a different subdirectory given by the `runName` paramter, in order to allow exploration of filtering parameters.

Filtering of counts matrix (to excluding noncoding transcripts and/or oscillating genes) was carried out with ***./filterCountsMatrix.R***

Running default DE analysis (excluding noncoding and oscilating genes) was carried out with ***./nf_differentialabundance.sh***

Alternatively, [RAPToR](https://github.com/LBMC/RAPToR/tree/master) can be run on the samples to determine developmental age and then DESeq2 is run manually with the ***./raptor_deseq2.R*** script.

## Custom plots

For both DEseq scripts above, the results tables are post processed to add genomic location and gene names using the ***./postProcessing*** script

Specific plots of gene expression on chrI and chrV are performed with the ***./plotResults.R***
