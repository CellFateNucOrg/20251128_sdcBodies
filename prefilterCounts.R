rm(list = ls())

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(yaml)

options(width=1e4)

serverPath="/Volumes/meister.data"
#serverPath="X:"

workDir=paste0(serverPath,"/sinem/sequencing/20261028_se_ill_rna_l3_wPM1556_autosomal_bmm1")
scriptDir=paste0(workDir,"/scripts")
setwd(workDir)

samplesheet<-read.csv(paste0(workDir,"/fileList.csv"), header=T, stringsAsFactors = F)


## params for filtering -----
outliers<-NULL  # vector of outlier samples or NULL


## gtf - load to find protein coding genes
genomeVer<-"WS298"
gtf<-import(paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))
gtf <- gtf[gtf$type == "gene" & gtf$gene_biotype=="protein_coding"]
mcols(gtf)<-mcols(gtf)[c("source","type","gene_id","gene_biotype","gene_name")]
gtf$source<-genomeVer
gtf<-sort(gtf)


## oscillating genes
#osc<-read.delim(paste0(serverPath,"/publicData/various/oscillatingGenes/oscillating_Hendriks_2014.tsv"))
#rising<-read.delim(paste0(serverPath,"/publicData/various/oscillatingGenes/rising_Hendriks_2014.tsv"))
osc1<-read.delim(paste0(serverPath,"/publicData/various/oscillatingGenes/oscillating_Meeuse_2020.tsv"))
osc2<-read.delim(paste0(serverPath,"/publicData/various/oscillatingGenes/oscillating_Latorre_2015.tsv"))

oscillating<-unique(c(osc1$wormbaseID,osc2$wormbaseID))


# filtering ------
## gene_counts
cnts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
dim(cnts)
if(!is.null(outliers)){
  idx<-colnames(cnts) %in% outliers
  noOut_cnts<-cnts[,!idx]
}

## gene_lengths
genelength<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_lengths.tsv"))
dim(genelength)

if(!is.null(outliers)){
  idx<-colnames(genelength) %in% outliers
  noOut_genelength<-genelength[,!idx]
}

## gene_tpm - needed for raptor (do not remove oscillating!)
genetpm<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_tpm.tsv"))
dim(genetpm)


## all genes but no outlier samples
if(!is.null(outliers)){
  idx<-colnames(cnts) %in% outliers
  noOut_cnts<-cnts[,!idx]
  write.table(noOut_cnts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_noOut.tsv"),sep="\t",row.names=F,quote=F)
  idx<-colnames(genelength) %in% outliers
  noOut_genelength<-genelength[,!idx]
  write.table(noOut_genelength,paste0(workDir,"/star_salmon/salmon.merged.gene_lengths_noOut.tsv"),sep="\t",row.names=F,quote=F)
  idx<-colnames(genetpm) %in% outliers
  noOut_genetpm<-genetpm[,!idx]
  write.table(noOut_genetpm,paste0(workDir,"/star_salmon/salmon.merged.gene_tpm_noOut.tsv"),sep="\t",row.names=F,quote=F)
  idx<-samplesheet$sample %in% outliers
  noOut_samplesheet<-samplesheet[!idx,]
  write.csv(noOut_samplesheet, paste0(workDir,"/samplesheet_noOut.csv"), row.names=F, quote=F)
}



## protein coding only
cnts<-cnts[cnts$gene_id %in% gtf$gene_id,]
genelength<-genelength[genelength$gene_id %in% gtf$gene_id,]
genetpm<-genetpm[genetpm$gene_id %in% gtf$gene_id,]
dim(cnts)
write.table(cnts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_pc.tsv"),sep="\t",row.names=F,quote=F)
write.table(genelength,paste0(workDir,"/star_salmon/salmon.merged.gene_lengths_pc.tsv"),sep="\t",row.names=F,quote=F)
write.table(genetpm,paste0(workDir,"/star_salmon/salmon.merged.gene_tpm_pc.tsv"),sep="\t",row.names=F,quote=F)

if(!is.null(outliers)){
  noOut_cnts<-noOut_cnts[noOut_cnts$gene_id %in% gtf$gene_id,]
  noOut_genelength<-noOut_genelength[noOut_genelength$gene_id %in% gtf$gene_id,]
  noOut_genetpm<-noOut_genetpm[noOut_genetpm$gene_id %in% gtf$gene_id,]
  dim(noOut_cnts)
  write.table(noOut_cnts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_pc_noOut.tsv"),sep="\t",row.names=F,quote=F)
  write.table(noOut_genelength,paste0(workDir,"/star_salmon/salmon.merged.gene_lengths_pc_noOut.tsv"),sep="\t",row.names=F,quote=F)
  write.table(noOut_genetpm,paste0(workDir,"/star_salmon/salmon.merged.gene_tpm_pc_noOut.tsv"),sep="\t",row.names=F,quote=F)
}


# remove oscillating
cnts<-cnts[!(cnts$gene_id %in% oscillating),]
genelength<-genelength[!(genelength$gene_id %in% oscillating),]
write.table(cnts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_pc_noOsc.tsv"),sep="\t",row.names=F,quote=F)
write.table(genelength,paste0(workDir,"/star_salmon/salmon.merged.gene_lengths_pc_noOsc.tsv"),sep="\t",row.names=F,quote=F)

if(!is.null(outliers)){
  noOut_cnts<-noOut_cnts[!(noOut_cnts$gene_id %in% oscillating),]
  noOut_genelength<-noOut_genelength[!(noOut_genelength$gene_id %in% oscillating),]
  dim(noOut_cnts)
  write.table(noOut_cnts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_pc_noOsc_noOut.tsv"),sep="\t",row.names=F,quote=F)
  write.table(noOut_genelength,paste0(workDir,"/star_salmon/salmon.merged.gene_lengths_pc_noOsc_noOut.tsv"),sep="\t",row.names=F,quote=F)
}



# create params.yaml --------
toCluster <- function(p) sub("^/Volumes/meister\\.data", "/mnt/meister.data", p)

minAbund=10L
minSamples=nrow(samplesheet)
suffixes <- c("", "_pc", "_pc_noOsc")

paramsets <- list(); i <- 0
for (suf in suffixes) for (shrink in c(FALSE, TRUE)) {
  i <- i + 1
  paramsets[[i]] <- list(
    paramset_name = sprintf("%02d_%s_minAbund%d_minSamples%d%s", i,
                            if (shrink) "lfcShrink" else "noShrink",
                            minAbund, minSamples, suf),
    deseq2_shrink_lfc = shrink,
    filtering_min_abundance = minAbund,
    filtering_min_samples = minSamples,
    matrix = file.path(toCluster(workDir), "star_salmon",
                       paste0("salmon.merged.gene_counts", suf, ".tsv")),
    transcript_length_matrix = file.path(toCluster(workDir), "star_salmon",
                                         paste0("salmon.merged.gene_lengths", suf, ".tsv"))
  )
}

write_yaml(paramsets, paste0(scriptDir,"/params.yml"),
           handlers = list(logical = verbatim_logical))


# create params_raptor.yaml --------
suffixes <- c("", "_pc")

paramsets_raptor <- list(); i <- 0
for (suf in suffixes) for (shrink in c(FALSE, TRUE)) {
  i <- i + 1
  paramsets_raptor[[i]] <- list(
    paramset_name = sprintf("%02d_%s_minAbund%d_minSamples%d%s_raptor", i,
                            if (shrink) "lfcShrink" else "noShrink",
                            minAbund, minSamples, suf),
    deseq2_shrink_lfc = shrink,
    filtering_min_abundance = minAbund,
    filtering_min_samples = minSamples,
    exploratory_main_variable = "group",
    matrix = file.path(toCluster(workDir), "star_salmon",
                       paste0("salmon.merged.gene_counts", suf, ".tsv")),
    transcript_length_matrix = file.path(toCluster(workDir), "star_salmon",
                                         paste0("salmon.merged.gene_lengths", suf, ".tsv"))
  )
}

write_yaml(paramsets_raptor, paste0(scriptDir,"/params_raptor.yml"),
           handlers = list(logical = verbatim_logical))


## param.yaml for noOutliers -----

if (length(outliers) > 0) {
  minAbund=10L
  minSamples=nrow(noOut_samplesheet)
  suffixes <- c("_noOut", "_pc_noOut", "_pc_noOsc_noOut")

  paramsets_noOut <- list(); i <- 0
  for (suf in suffixes) for (shrink in c(FALSE, TRUE)) {
    i <- i + 1
    paramsets_noOut[[i]] <- list(
      paramset_name = sprintf("%02d_%s_minAbund%d_minSamples%d%s", i,
                            if (shrink) "lfcShrink" else "noShrink",
                            minAbund, minSamples, suf),
      deseq2_shrink_lfc = shrink,
      filtering_min_abundance = minAbund,
      filtering_min_samples = minSamples,
      matrix = file.path(toCluster(workDir), "star_salmon",
                       paste0("salmon.merged.gene_counts", suf, ".tsv")),
      transcript_length_matrix = file.path(toCluster(workDir), "star_salmon",
                                         paste0("salmon.merged.gene_lengths",
                                                suf, ".tsv"))
    )
  }

  write_yaml(paramsets_noOut, paste0(scriptDir,"/params_noOut.yml"),
           handlers = list(logical = verbatim_logical))
}

## params_raptor_noOut.yml aptor runs without outliers -------
if (length(outliers) > 0) {
  suffixes <- c("_noOut", "_pc_noOut")

  paramsets_raptor_noOut <- list(); i <- 0
  for (suf in suffixes) for (shrink in c(FALSE, TRUE)) {
    i <- i + 1
    paramsets_raptor_noOut[[i]] <- list(
      paramset_name = sprintf("%02d_%s_minAbund%d_minSamples%d%s_raptor", i,
                              if (shrink) "lfcShrink" else "noShrink",
                              minAbund, minSamples, suf),
      deseq2_shrink_lfc = shrink,
      filtering_min_abundance = minAbund,
      filtering_min_samples = minSamples,
      exploratory_main_variable = "group",
      matrix = file.path(toCluster(workDir), "star_salmon",
                         paste0("salmon.merged.gene_counts", suf, ".tsv")),
      transcript_length_matrix = file.path(toCluster(workDir), "star_salmon",
                                           paste0("salmon.merged.gene_lengths",
                                                  suf, ".tsv"))
    )
  }

  write_yaml(paramsets_raptor_noOut, paste0(scriptDir,"/params_raptor_noOut.yml"),
             handlers = list(logical = verbatim_logical))
}



## names to pass to --paramset_name
paste(sapply(paramsets, `[[`, "paramset_name"), collapse = ",")
paste(sapply(paramsets_raptor, `[[`, "paramset_name"), collapse = ",")
if (length(outliers) > 0) {
  # names to pass to --paramset_name
  paste(sapply(paramsets_noOut, `[[`, "paramset_name"), collapse = ",")
  paste(sapply(paramsets_raptor_noOut, `[[`, "paramset_name"), collapse = ",")
}
