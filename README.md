# 20261103 autosomal_spots

## RNAseq alignment

Sample sheet and contrast matrix created with ***./makeSampleSheet.R***

Initial alignment to genome/transcriptome carried out with ***./nf_rnaseq_submit.sh***

## Differential expression analysis

Using params.yml to concurrently run multiple different standard configurations:
all data, only protein coding genes, removing oscillating genes +- lfcShrink and 
+- age estimation by [RAPToR](https://github.com/LBMC/RAPToR/tree/master) as a covariate (params_raptor.yml). It can also run analysis for batches where a particular outlier sample has been removed from analysis (params_noOut.yml).
The output is in subdirectories with the runName under each type of output. all in a single results folder (./diffabund)

Filtering of counts matrix (to excluding noncoding transcripts and/or oscillating genes) was carried out with ***./filterCountsMatrix.R***
This also prepares the params.yml

Then you must run ***./raptor_deseq2.R*** to get age estimates and add them to samplesheet_raptor.csv

You can then get a summary of the DESeq2 report by running da_summary/make_da_summary.sh. you need to specify which params*.yml files to use. 

## Custom plots

For both DEseq scripts above, the results tables are post processed to add genomic location and gene names using the ***./postProcessing*** script

Specific plots of gene expression on chrI and chrV are performed with the ***./plot_lfc_by_features_.R***
