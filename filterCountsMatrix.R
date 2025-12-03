library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)


source("./variables.R")
#genomeVer<-"WS295"

length(gtf)
table(gtf$type)
table(gtf$gene_biotype)
table(seqnames(gtf))

# protein coding only
#pc<-import(paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.protein_coding.gtf"))
if(samples=="no1298IPB2_lowInput"){
  counts<- read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.",samples,".tsv"))
} else {
  counts<- read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
}
dim(counts)
counts<-counts[counts$gene_id %in% gtf$gene_id,]
dim(counts)
#write.table(counts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_noRR_noSP.tsv"),sep="\t",row.names=F,quote=F)
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


# ## remove oscilating from combined bulk and low input sequencing
# counts<-read.delim(file="./star_salmon/salmon.merged.gene_counts.bulk_plus_lowinput.tsv",header=T)
#
# dim(counts)
# counts<-counts[counts$gene_id %in% pc$gene_id,]
# dim(counts)
#
# write.table(counts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts.bulk_plus_lowinput.protein_coding.tsv"),sep="\t",row.names=F,quote=F)
#
#
# counts_filt<-counts[!oscillating,]
# dim(counts_filt)
# dim(counts)
# write.table(counts_filt,paste0(workDir,"/star_salmon/salmon.merged.gene_counts.bulk_plus_lowinput.protein_coding_noOsc.tsv"),sep="\t",row.names=F,quote=F)
