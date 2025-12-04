# 20251128_sdcBodies

## RNAseq alignment

The easiest way to create a sample sheet is to first list the fastq files with full paths
and push them to a file which you then modify. The command will look something like this (depending on 
on the directory structure where your fastq files are stored and how they are called):

```
ls $PWD/fastq/*/*.fq.gz > fileList.txt
```

Then the sample sheet and contrast matrix are created with ***./makeSampleSheet.R***

NOTE: If you choose to filter out samples or add samples, you must make a new sample sheet, adding _<samples> to the name where <samples> is 
a name given to the sample set in variables.R (see below).

Initial alignment to genome/transcriptome carried out with ***./nf_rnaseq_submit.sh***

## Differential expression analysis

DESeq2 runs (with nf-core or raptor) is dependant on variables set in the ***./variables.R*** script.
For each combination of parameters, a unique results folder will be created.
Usually it is best to start with very permissive parameters and gradually perform stricter filtering if necessary.

If you choose to sensor particular samples, you must create a new sample sheet - add the code to ***./makeSampleSheet.R***
You do not need to filter columns of the counts matrix etc, as these will be selected according to the sample sheet, but the sample 
sheet should be named with fileList_<samples>.csv with samples being one of the variables set in the variables.R file.
If you combine with another dataset, you must add these samples not only to the sample sheet, but also to the
gene_counts, gene_lengths and gene_tpm files in star_salmon/. Add the code for doing that in ***./makeSampleSheet.R***, make sure to add
the _<samples> suffix to the name of those files too.

Results of DE analysis are each placed into a different subdirectory given by the `runName` paramter, in order to allow exploration of filtering parameters.
The each start with a prefix 'resNN' where NN is an ordinal number to keep track of the order of the analyses you perform.

Filtering of counts matrix (to exclude noncoding transcripts and/or oscillating genes) is carried out with ***./filterCountsMatrix.R***
Run this script even if you do not wish to filter any rows, so that the matrices will be named according to the input samples. (creates some
duplication but ensures you can go back to the relevant matrices even if you rerun nf_rnaseq_submit.sh which will overwrite the
default counts/lengths/tpm matrices).

Running default DE analysis (excluding noncoding and oscilating genes) is carried out with ***./nf_differentialabundance.sh***
If present, it uses ./variables.R to determine run parameters, otherwise you have to set them manually.

Alternatively, [RAPToR](https://github.com/LBMC/RAPToR/tree/master) can be run on the samples to determine developmental age and then DESeq2 is run on the data using the age as a covariate with the ***./raptor_deseq2.R*** script.

## Custom plots

These scripts also depend on variables in the variables.R file. Make sure you have the right combination of parameters entered to the scripts choose the correct run directory.

For both DEseq scripts above, the results tables are post processed to add genomic location and gene names using the ***./postProcessing*** script

Specific plots of gene expression on chrI and chrV are performed with the ***./plotResults.R***
