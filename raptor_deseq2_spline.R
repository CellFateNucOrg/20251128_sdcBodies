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
source("./filterCountsMatrix.R")

print(runName)

source(paste0(workDir,"/raptor_functions.R"))
col.palette1 <- c('grey20', 'firebrick', 'royalblue', 'forestgreen')


contrasts<-read.csv(contrastsFile,sep=",",header=T)
counts<-read.delim(countsFile)
tpm<-read.delim(tpmFile)
sampleSheet<-read.csv(sampleSheetFile, header=T, stringsAsFactors = T)

dir.create(paste0(workDir,runName,"/plots"), recursive=T, showWarnings = FALSE)
dir.create(paste0(workDir,runName,"/tables/differential"), recursive=T, showWarnings = FALSE)


# load reference
ref <- RAPToR::prepare_refdata("Cel_larval", "wormRef", 600)

qs=list()

counts<-counts[rowSums(counts[, c(-1,-2)]>minAbund)>minSamples,]

qs$count<-round(counts[, c(-1,-2)],0)
row.names(qs$count)<-counts$gene_id
colnames(qs$count)<-gsub("^X","",colnames(qs$count)) # replace X at beginning of names that start with numbers

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

ae_X <- RAPToR::ae(qs$tpm, ref,prior=25, prior.params=3)

print(ae_X)
par(mfrow=c(1,1))

# plot ages
lastSample=nrow(ae_X$age.estimates)
pdf(paste0(workDir,runName,"/plots/",prefix,"ages.pdf"),width=9,height=6*(floor(lastSample/20)+1))
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
lastSample=nrow(ae_X$age.estimates)
for (i in seq(1,lastSample,by=6)) {
  plot_cor(ae_X,subset=i:min(i+5,lastSample))
}
dev.off()
par(mfrow=c(1,1))


qs_rc <- RAPToR::ref_compare(
  X = log1p(qs$tpm), # sample data, log(TPM+1)
  ref = ref, # ref object
  ae_obj = ae_X, # ae object
  group = qs$pdat$group # factor defining compared groups (wt/mut)
)


print(qs_rc)
plotList=list()
c=1
for(c in 1:nrow(contrasts)){
  qs_lfc<-RAPToR::get_logFC(qs_rc,l=contrasts$target[c],l0=contrasts$reference[c],
                    verbose=T)

  plotList[[contrasts$id[c]]]<-gg_logFC(qs_lfc, main = contrasts$id[c],
               xlab = contrasts$id[c], ylab = "Developmental log2FC")
}


p<-ggarrange(plotlist=plotList,ncol=3)
p<-annotate_figure(p, top=text_grob(runName))
ggsave(paste0(workDir,runName,"/plots/",prefix,"_development_lfc.pdf"),p,height=10, width=29, units="cm")

ageRange<-range(ae_X$age.estimates[,1] + c(-1,1))
dfs<-2:8

dfs_ssq<-find_df(ref,w=ageRange,dfs=dfs)

plot(dfs, dfs_ssq, type='b', xlab="Degrees of freedom", ylab="SSQ per sample")

df_opti<-5

g.filt <- RAPToR::format_to_ref(qs$count, refdata = ref$interpGE)$samp

idx<-match(qs$pdat$sample,row.names(ae_X$age.estimates))
qs$pdat$age<-ae_X$age.estimates[idx,1]
qs$pdat$ageZnorm<- (qs$pdat$age-mean(qs$pdat$age))/sd(qs$pdat$age)

DEref.shifts <- run_DESeq2_ref(X = g.filt, # sample counts
               p = qs$pdat, # sample pdata
               ae_values = qs$pdat$age, # age estimates
               formula = ~group, # formula (age is added by the function)
               ref = ref, # ref object
               ns.df = df_opti, # df for spline fit
               window.extend = 1) # reference window extension


c=3
for(c in 1:nrow(contrasts)){
  res <- results(DEref.shifts,contrast=c("group",contrasts$target[c],
                                         contrasts$reference[c]),alpha=0.05)
  if(shrink){
    res <- lfcShrink(DEref.shifts,contrast=c("group",contrasts$target[c],
                                             contrasts$reference[c]),alpha=0.05,type="ashr")
  }
  summary(res)
  write.table(res, file=paste0(workDir,runName,"/tables/differential/",prefix,contrasts$id[c],".deseq2.results.tsv"), sep="\t", quote=F)
}

# # Get results table from deseq output, managing NAs
# if(!is.null(coefname)){
#   res <- results(dds, name=coefname,alpha=0.05)
# } else {
#   res <- results(dds, contrast=contrastname, alpha=0.05)
# }
# # manage NAs
# res$padj[is.na(res$padj)] <- 1
# res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
# return(res)
#
#
#
#
#
#
#
#
# ae_X
#
# qs$count <- qs$count[apply(qs$count, 1, max)>5, ]
#
# # use age instead of batch
# dd0 <- DESeqDataSetFromMatrix(countData = qs$count,
#                               colData = qs$pdat,
#                               design = ~ageZnorm+group)
# dd0<-DESeq(dd0)
#
#
# for(c in 1:nrow(contrasts)){
#   res <- results(dd0,contrast=c("group",contrasts$target[c],contrasts$reference[c]),alpha=0.05)
#   if(shrink){
#   res <- lfcShrink(dd0,contrast=c("group",contrasts$target[c],contrasts$reference[c]),alpha=0.05,type="ashr")
#   }
#   summary(res)
#   write.table(res, file=paste0(workDir,runName,"/tables/differential/",prefix,contrasts$id[c],".deseq2.results.tsv"), sep="\t", quote=F)
# }
#
#
#
#
#
# vsd <- vst(dd0, blind=FALSE)
# p1<-plotPCA(vsd, intgroup=c("group"))
# p2<-plotPCA(vsd, intgroup=c("group"),pcsToUse=3:4)
# p3<-plotPCA(vsd, intgroup=c("group"),pcsToUse=5:6)
# p<-ggarrange(p1,p2,p3,nrow=3)
# ggsave(paste0(workDir,runName,"/plots/",prefix,"_PCA.pdf"),p,height=29, width=19, units="cm")
#

file.copy(from=paste0(workDir,"/raptor_deseq2_spline.R"),
          to=paste0(workDir,runName,"/raptor_deseq2_spline.R"),
          overwrite=TRUE)

file.copy(from=paste0(workDir,"/variables.R"),
          to=paste0(workDir,runName,"/variables.R"),
          overwrite=TRUE)

devtools::session_info(to_file=paste0(workDir,runName,"/logs/sessionInfo_raptor_deseq2.txt"))

