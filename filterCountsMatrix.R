library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)


source("./variables.R")

gtf<-import(gtfFile)
length(gtf)
table(gtf$type)
table(gtf$gene_biotype)
table(seqnames(gtf))

# import the right unfiltered count matrix depending on samples used
if(samples=="allSamples"){
  counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
} else if(samples=="no1298IPB2"){
  #print("no1298IPB2")
  counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
} else if(samples=="no1298IPB2_lowInput"){
  #print("no1298IPB2_lowInput")
  counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.no1298IPB2_lowInput.tsv"))
}


dim(counts)
if(filterNoncoding){
  counts<-counts[counts$gene_id %in% gtf$gene_id,]
  dim(counts)
}

if(filterMitochondrial){
  idx<-seqnames(gtf)=="MtDNA"
  if(sum(idx)>0){
    mitoGenes<-gtf$gene_id[idx]
    counts<-counts[!counts$gene_id %in% mitoGenes,]
  }
  dim(counts)
}
write.table(counts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts",
                          ifelse(filterNoncoding,".protein_coding",""),
                          ifelse(filterMitochondrial,".noMito",""),
                          ifelse(samples=="allSamples","",paste0(".",samples)),
                          ".tsv"),sep="\t",row.names=F,quote=F)

# remove oscillating genes
latorre<-read.delim(paste0(serverPath,"/publicData/Various/Oscillating_genes/oscillatingGenes_latorre.tsv"))

meeuse<-read.delim(paste0(serverPath,"/publicData/Various/Oscillating_genes/oscillatingGenes_Meeuse.tsv"))

oscillating<-counts$gene_id %in% latorre$wormbaseID | counts$gene_id %in% meeuse$wormbaseID
counts_filt<-counts[!oscillating,]
dim(counts_filt)
dim(counts)
write.table(counts_filt,paste0(workDir,"/star_salmon/salmon.merged.gene_counts",
                          ifelse(filterNoncoding,".protein_coding",""),
                          ifelse(filterMitochondrial,".noMito",""),
                          ifelse(filterOscillating,".noOsc",""),
                          ifelse(samples=="allSamples","",paste0(".",samples)),
                          ".tsv"),sep="\t",row.names=F,quote=F)

