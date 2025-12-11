library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)


source("./variables.R")


gtf<-import(gtfFile)
gtf <- gtf[gtf$type == "gene"]
length(gtf)

#####################-
# filter counts matrix -----
#####################-
# import the right unfiltered count matrix depending on samples used
if(samples=="allSamples"){
  counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
# } else if(samples=="no1298IPB2"){
#   #print("no1298IPB2")
#   counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
# } else if(samples=="no1298IPB2_lowInput"){
#   #print("no1298IPB2_lowInput")
#   counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.no1298IPB2_lowInput.tsv"))
} else if(samples=="allSamples_lowInput"){
  #print("no1298IPB2_lowInput")
  counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.allSamples_lowInput.tsv"))
}

#left_join(counts,data.frame(gtf),by=c("gene_id"="gene_id"))


# remove noncoding genes
dim(counts)
if(filterNoncoding){
  counts<-counts[counts$gene_id %in% gtf$gene_id,]
  dim(counts)
}

# remove mitochondrial genes
if(filterMitochondrial){
  idx<-seqnames(gtf)=="MtDNA"
  if(sum(idx)>0){
    mitoGenes<-gtf$gene_id[idx]
    counts<-counts[!counts$gene_id %in% mitoGenes,]
  }
  dim(counts)
}

# remove oscillating genes
if(filterOscillating){

  latorre<-read.delim(paste0(serverPath,"/publicData/Various/Oscillating_genes/oscillatingGenes_latorre.tsv"))
  meeuse<-read.delim(paste0(serverPath,"/publicData/Various/Oscillating_genes/oscillatingGenes_Meeuse.tsv"))

  oscillating<-counts$gene_id %in% latorre$wormbaseID | counts$gene_id %in% meeuse$wormbaseID
  counts<-counts[!oscillating,]
  dim(counts)
  dim(counts)
}

outFile=paste0(workDir,"/star_salmon/salmon.merged.gene_counts",
                   ifelse(filterNoncoding,".protein_coding",""),
                   ifelse(filterMitochondrial,".noMito",""),
                   ifelse(filterOscillating,".noOsc",""),
                   ifelse(samples=="allSamples","",paste0(".",samples)),
                   ".tsv")

# never overwrite the original file
if (outFile!=paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv")){
  print(paste0("writing ",outFile))
  write.table(counts,outFile,sep="\t",row.names=F,quote=F)
}



