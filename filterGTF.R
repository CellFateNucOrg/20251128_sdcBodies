library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

source("./variables.R")

Celegans_WB<-Celegans
seqlevelsStyle(Celegans_WB)<-"NCBI"
seqlevels(Celegans_WB)<-c("I", "II", "III", "IV", "V", "X", "MtDNA")


####################-
# filter gtf ------
###################-

gtf<-import(paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))
length(gtf)
table(gtf$type)
table(gtf$gene_biotype)
table(seqnames(gtf))
gtf<-gtf[gtf$type=="gene"]
length(gtf)

seqlevels(gtf)<-c("I", "II", "III", "IV", "V", "X", "MtDNA")
seqlengths(gtf)<-seqlengths(Celegans_WB)
gtf<-sort(gtf)
table(seqnames(gtf))
print(length(gtf))


# filter for protein coding genes
filterNoncoding=T
filterMitochondrial=F
idx<-gtf$gene_biotype=="protein_coding"
gtf_filt<-gtf[idx,]
table(gtf_filt$gene_biotype)
print(length(gtf_filt))

# export filtered gtf
rtracklayer::export(gtf_filt,con=paste0(serverPath,"/publicData/genomes/",genomeVer,
                                   "/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset",
                                   ifelse(filterNoncoding,".protein_coding",""),
                                   ifelse(filterMitochondrial,".noMito",""),
                                   ".gtf"))

# filter for protein coding and nuclear genes
filterNoncoding=T
filterMitochondrial=T
idx<-gtf$gene_biotype=="protein_coding"
gtf_filt<-gtf[idx,]
idx<-seqnames(gtf_filt)!="MtDNA"
gtf_filt<-gtf_filt[idx,]
table(seqnames(gtf_filt))
print(length(gtf_filt))

# export filtered gtf
rtracklayer::export(gtf_filt,con=paste0(serverPath,"/publicData/genomes/",genomeVer,
                                        "/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset",
                                        ifelse(filterNoncoding,".protein_coding",""),
                                        ifelse(filterMitochondrial,".noMito",""),
                                        ".gtf"))






