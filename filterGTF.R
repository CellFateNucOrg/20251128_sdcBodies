library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(GenomeInfoDb)

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

seqlevels(gtf)<-c("I", "II", "III", "IV", "V", "X", "MtDNA")
seqlengths(gtf)<-seqlengths(Celegans_WB)
gtf<-sort(gtf)
table(seqnames(gtf))
print(length(gtf))

gtf_gene<-gtf[gtf$type=="gene"]
print(length(gtf_gene))




# filter for protein coding genes
filterNoncoding=T
filterMitochondrial=F
idx<-gtf_gene$gene_biotype=="protein_coding"
gtf_gene_filt<-gtf_gene[idx,]
table(gtf_gene_filt$gene_biotype)
print(length(gtf_gene_filt))

print(length(gtf[gtf$gene_id %in% gtf_gene_filt$gene_id]))

# export filtered gtf
rtracklayer::export(gtf[gtf$gene_id %in% gtf_gene_filt$gene_id],
                    con=paste0(serverPath,"/publicData/genomes/",genomeVer,
                               "/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset",
                               ifelse(filterNoncoding,".protein_coding",""),
                               ifelse(filterMitochondrial,".noMito",""),
                               ".gtf"))

# filter for protein coding and nuclear genes
filterNoncoding=T
filterMitochondrial=T
idx<-gtf_gene$gene_biotype=="protein_coding"
gtf_gene_filt<-gtf_gene[idx,]
idx<-seqnames(gtf_gene_filt)!="MtDNA"
gtf_gene_filt<-gtf_gene_filt[idx,]
table(seqnames(gtf_gene_filt))
print(length(gtf_gene_filt))

print(length(gtf[gtf$gene_id %in% gtf_gene_filt$gene_id]))

# export filtered gtf
rtracklayer::export(gtf[gtf$gene_id %in% gtf_gene_filt$gene_id],
                    con=paste0(serverPath,"/publicData/genomes/",genomeVer,
                               "/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset",
                               ifelse(filterNoncoding,".protein_coding",""),
                               ifelse(filterMitochondrial,".noMito",""),
                               ".gtf"))



