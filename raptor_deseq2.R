library(RAPToR)
library(wormRef)
library(rtracklayer)
library(DESeq2)
library(splines)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(plotly)
library(ggrepel)
library(rstatix)


theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)

serverPath="/Volumes/meister.data"
#serverPath="Z:/MeisterLab"
workDir=paste0(serverPath,"/kalyan/sequencing/20261028_kg_ill_rna_l3_wPM1486_autosomal_smc1")

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
samplesheet <- read.csv(file.path(workDir, "samplesheet.csv"),header=T)
minAbund=10
minSamples=nrow(samplesheet)
outliers<-c("IPTG_1486_3_L2","IPTG_1486_3_L7") #NULL  # vector of outlier samples or NULL

# age estimation parameters:
prior=17.5 #NULL # expected age in hours (at 20C) - try 27-28 since we grow at 22C, or look at plots with NULL first
prior.params=10 #NULL # wide SD, usually i use 10 to help confine it to same larval stage, but not bias it too much. use NULL for initial plot
refName<-"Cel_larval" # could also try "Cel_larv_YA" - seems better age estimation with NULL prior, but worse corelation overall.

runName=paste0("01_noShrink_minAbund",minAbund,"_minSamples",minSamples)

prefix=""
genomeVer<-"WS298"

setwd(workDir)

source(paste0(workDir,"/scripts/raptor_functions.R"))
col.palette1 <- c('grey20', 'firebrick', 'royalblue', 'forestgreen')
#mapMethod="salmon"
#prefix=paste0(mapMethod,"_raptor_age_deseq2_pc_nomito_lfcShrink_")



counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
tpm<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_tpm.tsv"))


ss<-read.csv(paste0(workDir,"/samplesheet.csv"), header=T, stringsAsFactors = T)


idx<-which(colnames(counts) %in% ss$sample)
counts<-counts[,c(1,2,idx)]
tpm<-tpm[,c(1,2,idx)]



shrink=ifelse(grepl("noShrink",runName),FALSE,TRUE)
if(shrink==TRUE){
  ylimits=c(-0.5,0.5)
} else {
  ylimits=c(-2,2)
}



dir.create(paste0(workDir,"/diffabund/custom/raptor/",runName), recursive=T, showWarnings = FALSE)

#dir.create(paste0(workDir,"/diffabund/custom/differential"), recursive=T, showWarnings = FALSE)

gtf<-import(paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.protein_coding.gtf"))

plot_refs("wormRef")

# load reference
ref <- prepare_refdata(refName, "wormRef", 600)


#pdat<-ss[,c("sample","genotype","strain","group","replicate","group")]

qs=list()

counts<-counts[rowSums(counts[, c(-1,-2)]>=minAbund)>=minSamples,]

qs$count<-round(counts[, c(-1,-2)],0)
row.names(qs$count)<-counts$gene_id
colnames(qs$count)<-gsub("^X","",colnames(qs$count))

#tpm<-read.delim(paste0("./",ifelse(mapMethod=="star","star_",""),"salmon/salmon.merged.gene_tpm.tsv"),header=T)
tpm<-tpm[tpm$gene_id %in% row.names(qs$count),]
X<- tpm[, c(-1,-2)]
row.names(X)<-tpm$gene_id
colnames(X)<-gsub("^X","",colnames(X))
qs$tpm<-X

pdatOrder<-match(colnames(X),ss$sample)
qs$pdat<-ss[pdatOrder,]

str(qs$pdat)
qs$pdat$group<-factor(qs$pdat$group)
qs$pdat$strain<-factor(qs$pdat$strain)
qs$pdat$replicate<-factor(qs$pdat$replicate)

summary(qs$pdat)

ae_X <- ae(qs$tpm, ref, prior=prior, prior.params=prior.params)
print(ae_X)

par(mfrow=c(1,1))

# plot ages
pdf(paste0(workDir,"/diffabund/custom/raptor/",runName,"/",refName,"_prior",prior,"_ages.pdf"),width=9,height=6)
plot(ae_X,group=qs$pdat$strain, main="Ages of samples by strain",
     xlab="Age (h)", cex=0.5, col=col.palette1[qs$pdat$replicate])
plot(ae_X,group=qs$pdat$group,main="Ages of samples by group",
     xlab="Age (h)", cex=0.5, col=col.palette1[qs$pdat$replicate])
plot(ae_X,group=qs$pdat$replicate, main="Ages of samples by replicate",
     xlab="Age (h)", cex=0.5, col=col.palette1[qs$pdat$group])
dev.off()

# plot full correlation
pdf(paste0(workDir,"/diffabund/custom/raptor/",runName,"/",refName,"_prior",prior,"_full_age_correlations.pdf"),width=9,height=6)
par(mfrow=c(2,3))
pages <- split(seq_len(nrow(samplesheet)), ceiling(seq_len(nrow(samplesheet)) / 6))
for (s in pages) plot_cor(ae_X, subset = s)
dev.off()
par(mfrow=c(1,1))


qs_rc <- ref_compare(
  X = log1p(qs$tpm), # sample data, log(TPM+1)
  ref = ref, # ref object
  ae_obj = ae_X, # ae object
  group = qs$pdat$group # factor defining compared groups (wt/mut)
)


print(qs_rc)
plotList=list()
c=1
for(c in 1:nrow(contrasts)){
  qs_lfc<-get_logFC(qs_rc,l=contrasts$target[c],l0=contrasts$reference[c],
                    verbose=T)

  plotList[[contrasts$id[c]]]<-gg_logFC(qs_lfc, main = contrasts$id[c],
               xlab = contrasts$id[c], ylab = "Developmental log2FC")
}


p<-ggarrange(plotlist=plotList,ncol=nrow(contrasts))
ggsave(paste0(workDir,"/diffabund/custom/raptor/",runName,"/",refName,"_prior",prior,"_development_lfc.pdf"),p,height=10, width=29, units="cm")



idx<-match(qs$pdat$sample,row.names(ae_X$age.estimates))
qs$pdat$age<-ae_X$age.estimates[idx,1]
qs$pdat$ageZnorm<- (qs$pdat$age-mean(qs$pdat$age))/sd(qs$pdat$age)

write.csv(qs$pdat,paste0(workDir,"/samplesheet_raptor.csv"),
          row.names=F, quote=F)

if(!is.null(outliers)){
  idx<-qs$pdat$sample %in% outliers
  noOutss<-qs$pdat[!idx,]
  write.csv(noOutss,paste0(workDir,"/samplesheet_raptor_noOut.csv"),
            row.names=F, quote=F)
}

##### optional ######


ae_X

qs$count <- qs$count[apply(qs$count, 1, max)>5, ]
# use age instead of batch
dd0 <- DESeqDataSetFromMatrix(countData = qs$count,
                              colData = qs$pdat,
                              design = ~ageZnorm+group)
dd0<-DESeq(dd0)


for(c in 1:nrow(contrasts)){
  res <- results(dd0,contrast=c("group",contrasts$target[c],contrasts$reference[c]),alpha=0.05)
  if(shrink){
    res <- lfcShrink(dd0,contrast=c("group",contrasts$target[c],contrasts$reference[c]),alpha=0.05,type="ashr")
  }
  summary(res)
  write.table(res, file=paste0(workDir,"/diffabund/custom/raptor/",runName,"/",contrasts$id[c],".deseq2.results.tsv"), sep="\t", quote=F)
}


vsd <- vst(dd0, blind=FALSE)
p1<-plotPCA(vsd, intgroup=c("group"))
p2<-plotPCA(vsd, intgroup=c("group"),pcsToUse=3:4)
p3<-plotPCA(vsd, intgroup=c("group"),pcsToUse=5:6)
p<-ggarrange(p1,p2,p3,nrow=3)
ggsave(paste0(workDir,"/diffabund/custom/raptor/PCA.pdf"),p,height=29, width=19, units="cm")


file.copy(from=paste0(workDir,"/raptor_deseq2.R"),
          to=paste0(workDir,"/diffabund/custom/raptor/raptor_deseq2.R"),
          overwrite=TRUE)
