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

source("./variables.R")

#serverPath="/Volumes/external.data/MeisterLab"
#serverPath="Z:/MeisterLab"
#workDir=paste0(serverPath,"/jsemple/20251118_sdcBodies")
#minAbund=10
#minSamples=16
#runName=paste0("/res11_no1298IPB2_minAbund",minAbund,"_minSamples",minSamples,"_pc_raptor")
#runName=paste0("/res12_no1298IPB2_minAbund",minAbund,"_minSamples",minSamples,"_pc_raptor")
#runName=paste0("/res13_no1298IPB2_minAbund",minAbund,"_minSamples",minSamples,"_pc_raptor_noShrink")
#runName=paste0("/res14_no1298IPB2_minAbund",minAbund,"_minSamples",minSamples,"_pc_raptor_noShrink")
#runName=paste0("/res19_no1298IPB2_lowInput_minAbund",minAbund,"_minSamples",minSamples,"_pc_raptor_noShrink")
#runName=paste0("/res20_no1298IPB2_minAbund",minAbund,"_minSamples",minSamples,"_pc_raptor_noShrink")
#prefix=""
#genomeVer<-"WS295"

#setwd(workDir)

source(paste0(workDir,"/raptor_functions.R"))
col.palette1 <- c('grey20', 'firebrick', 'royalblue', 'forestgreen')
#mapMethod="salmon"
#prefix=paste0(mapMethod,"_raptor_age_deseq2_pc_nomito_lfcShrink_")

#contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)

#counts<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.protein_coding.tsv"))
#tpm<-read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_tpm.tsv"))


#sampleSheet<-read.csv(paste0(workDir,"/fileList_extended_no1298IPB2.csv"), header=T, stringsAsFactors = T)


idx<-which(colnames(counts) %in% sampleSheet$sample)
counts<-counts[,c(1,2,idx)]
tpm<-tpm[,c(1,2,idx)]



shrink=ifelse(grepl("noShrink",runName),FALSE,TRUE)
if(shrink==TRUE){
  ylimits=c(-0.5,0.5)
} else {
  ylimits=c(-2,2)
}



dir.create(paste0(workDir,runName,"/plots"), recursive=T, showWarnings = FALSE)

dir.create(paste0(workDir,runName,"/tables/differential"), recursive=T, showWarnings = FALSE)

#gtf<-import(paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.protein_coding.gtf"))


# load reference
ref <- prepare_refdata("Cel_larval", "wormRef", 600)

#pdat<-sampleSheet[,c("sample","genotype","strain","group","replicate","group")]

qs=list()

counts<-counts[rowSums(counts[, c(-1,-2)]>minAbund)>minSamples,]

qs$count<-round(counts[, c(-1,-2)],0)
row.names(qs$count)<-counts$gene_id
colnames(qs$count)<-gsub("^X","",colnames(qs$count))

#tpm<-read.delim(paste0("./",ifelse(mapMethod=="star","star_",""),"salmon/salmon.merged.gene_tpm.tsv"),header=T)
tpm<-tpm[tpm$gene_id %in% row.names(qs$count),]
X<- tpm[, c(-1,-2)]
row.names(X)<-tpm$gene_id
colnames(X)<-gsub("^X","",colnames(X))
qs$tpm<-X

pdatOrder<-match(colnames(X),sampleSheet$sample)
qs$pdat<-sampleSheet[pdatOrder,]

str(qs$pdat)
qs$pdat$group<-factor(qs$pdat$group)
qs$pdat$strain<-factor(qs$pdat$strain)
qs$pdat$replicate<-factor(qs$pdat$replicate)

summary(qs$pdat)

ae_X <- ae(qs$tpm, ref)

print(ae_X)
par(mfrow=c(1,1))

# plot ages
pdf(paste0(workDir,runName,"/plots/",prefix,"ages.pdf"),width=9,height=6)
plot(ae_X,group=qs$pdat$strain, main="Ages of samples by strain",
     xlab="Age (h)", cex=0.5, col=col.palette1[qs$pdat$replicate])
plot(ae_X,group=qs$pdat$group,main="Ages of samples by group",
     xlab="Age (h)", cex=0.5, col=col.palette1[qs$pdat$replicate])
plot(ae_X,group=qs$pdat$replicate, main="Ages of samples by replicate",
     xlab="Age (h)", cex=0.5, col=col.palette1[qs$pdat$group])
dev.off()

# plot full correlation
pdf(paste0(workDir,runName,"/plots/",prefix,"_full_age_correlations.pdf"),width=9,height=6)
par(mfrow=c(2,3))
plot_cor(ae_X,subset=1:6)
plot_cor(ae_X,subset=7:12)
plot_cor(ae_X,subset=13:17)
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


p<-ggarrange(plotlist=plotList,ncol=3)
ggsave(paste0(workDir,runName,"/plots/",prefix,"_development_lfc.pdf"),p,height=10, width=29, units="cm")



idx<-match(qs$pdat$sample,row.names(ae_X$age.estimates))
qs$pdat$age<-ae_X$age.estimates[idx,1]
qs$pdat$ageZnorm<- (qs$pdat$age-mean(qs$pdat$age))/sd(qs$pdat$age)

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
  write.table(res, file=paste0(workDir,runName,"/tables/differential/",prefix,contrasts$id[c],".deseq2.results.tsv"), sep="\t", quote=F)
}





vsd <- vst(dd0, blind=FALSE)
p1<-plotPCA(vsd, intgroup=c("group"))
p2<-plotPCA(vsd, intgroup=c("group"),pcsToUse=3:4)
p3<-plotPCA(vsd, intgroup=c("group"),pcsToUse=5:6)
p<-ggarrange(p1,p2,p3,nrow=3)
ggsave(paste0(workDir,runName,"/plots/",prefix,"_PCA.pdf"),p,height=29, width=19, units="cm")


file.copy(from=paste0(workDir,"/raptor_deseq2.R"),
          to=paste0(workDir,runName,"/raptor_deseq2.R"),
          overwrite=TRUE)

file.copy(from=paste0(workDir,"/variables.R"),
          to=paste0(workDir,runName,"/variables.R"),
          overwrite=TRUE)

devtools::session_info(to_file=paste0(workDir,runName,"/logs/sessionInfo_raptor_deseq2.txt"))
