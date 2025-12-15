library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(plotly)
library(ggrepel)
library(rstatix)
library(GenomeInfoDb)
library(plotgardener)
library(TxDb.Celegans.UCSC.ce11.refGene)
library(org.Ce.eg.db)

theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)


#source("./variables.R")
runName<-"/res05_allSamples_lowInput_minAbund10_minSamples17_pc_noMito_raptor_noShrink"
serverPath="/Volumes/external.data/MeisterLab"
workDir=paste0(serverPath,"/jsemple/20251118_sdcBodies")
source("./functions_exploratory_analysis.R")
prefix=""

contrasts<-read.csv(contrastsFile,sep=",",header=T)


results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir, "/forPaper"), showWarnings = FALSE, recursive = TRUE)


### binned lfc ------
gr<-tableToGranges(results)
seqlevelsStyle(gr)<-"UCSC"

targets<-import(paste0(workDir,"/tracks/targets.bed"))

binWidth=3e5
tiles<-tileGenome(seqlengths(Celegans),tilewidth=binWidth, cut.last.tile.in.chrom = T)
tiles$target<-F
ol<-findOverlaps(tiles,resize(targets,width=1,fix="center"))
tiles$target[queryHits(ol)]<-T

#tiles$binNum<-1:length(tiles)
ol<-findOverlaps(resize(gr,width=1,fix="center"),tiles)
gr$binNum<-subjectHits(ol)
gr$target<-tiles$target[subjectHits(ol)]

df<-data.frame(gr)
df$binNum<-factor(df$binNum)

ylim=c(-1.3,1.3)

# chrI
minBin=min(as.numeric(df$binNum[df$seqnames=="chrI"]))
maxBin=max(as.numeric(df$binNum[df$seqnames=="chrI"]))

subdf<-df[df$seqnames=="chrI" & df$group=="PMW1001_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)


p1<-ggplot(subdf, aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target),outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("PMW1001 [ChrV-lacO, GBPD1-lacI::mCherry, SDC-1::GFP]")) +
  geom_hline(yintercept=0,alpha=0.3) +
  scale_x_discrete(name="ChrI (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(minBin,maxBin,5)*binWidth/1e6) +
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90, size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3) +
  ylab("Log<sub>2</sub>FC ( Bound / IPTG )") +
  theme(plot.title=element_text(face="bold"))
p1



subdf<-df[df$seqnames=="chrI" & df$group=="PMW1291_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

p1a<-ggplot(subdf,aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target),outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("PMW1291 [ChrI-lacO, GBPD1-lacI::mCherry, SDC-1::GFP]")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrI (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(minBin,maxBin,5)*binWidth/1e6) +
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90,size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)+
  ylab("Log<sub>2</sub>FC ( Bound / IPTG )")+
  theme(plot.title=element_text(face="bold"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p1a



subdf<-df[df$seqnames=="chrI" & df$group=="PMW1298_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

p1b<-ggplot(subdf,aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target),outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("PMW1298 [ChrI-lacO, lacI::mCherry, SDC-1::GFP]")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrI genomic position (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(minBin,maxBin,5)*binWidth/1e6) +
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90, size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)+
  ylab("Log<sub>2</sub>FC ( Bound / IPTG )")+
  theme(plot.title=element_text(face="bold"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p1b

p<-ggarrange(p1a,p1b,p1,nrow=3)
p
p<-annotate_figure(p, top=text_grob("ChrI binned LFC",face="bold",size=14))
ggsave(paste0(workDir,"/forPaper/","Fig_chrI_",binWidth/1e3,"kb.pdf"),p,width=27,height=19,unit="cm")


# chrV

minBin=min(as.numeric(df$binNum[df$seqnames=="chrV"]))
maxBin=max(as.numeric(df$binNum[df$seqnames=="chrV"]))

subdf<-df[df$seqnames=="chrV" & df$group=="PMW1001_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

p2<-ggplot(subdf, aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target),outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("PMW1001 [ChrV-lacO, GBPD1-lacI::mCherry, SDC-1::GFP]")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrV genomic position (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(1,maxBin-minBin-1,5)*binWidth/1e6)+
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90,size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)+
  ylab("Log<sub>2</sub>FC ( Bound / IPTG )")+
  theme(plot.title=element_text(face="bold"))
p2



subdf<-df[df$seqnames=="chrV" & df$group=="PMW1291_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

p2a<-ggplot(subdf, aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target), outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("PMW1291 [ChrI-lacO, GBPD1-lacI::mCherry, SDC-1::GFP]")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrV genomic position (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(1,maxBin-minBin-1,5)*binWidth/1e6)+
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90, size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)+
  ylab("Log<sub>2</sub>FC ( Bound / IPTG )") +
  theme(plot.title=element_text(face="bold"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2a



subdf<-df[df$seqnames=="chrV" & df$group=="PMW1298_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

p2b<-ggplot(subdf,aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target), outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("PMW1298 [ChrI-lacO, lacI::mCherry, SDC-1::GFP]")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrV genomic position (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(1,maxBin-minBin-1,5)*binWidth/1e6)+
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90, size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)+
  ylab("Log<sub>2</sub>FC ( Bound / IPTG )") +
  theme(plot.title=element_text(face="bold"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2b

p<-ggarrange(p2a,p2b,p2,nrow=3)
p
p<-annotate_figure(p, top=text_grob("ChrV binned LFC",face="bold",size=14))
ggsave(paste0(workDir,"/forPaper/","Fig_chrV_",binWidth/1e3,"kb.pdf"),p,width=27,height=19,unit="cm")





###### zoomed lineplots -----
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
gr<-tableToGranges(results)
seqlevelsStyle(gr)<-"UCSC"

targets<-import(paste0(workDir,"/tracks/targets.bed"))
seqinfo(targets)
seqlevels(targets)<-seqlevels(Celegans)
seqinfo(targets)<-seqinfo(Celegans)

regionSize=3e5

targetRegions<-trim(resize(targets,width=regionSize,fix="center"))

df<-data.frame(subsetByOverlaps(gr,targetRegions))
targetdf<-data.frame(targets)
targetdf$log2FoldChange<-0

ylim=c(-1.5,0.6)

df$group<-factor(df$group, levels=c("PMW1291_Bound_vs_IPTG",
                             "PMW1298_Bound_vs_IPTG",
                             "PMW1001_Bound_vs_IPTG"),
                 labels=c("PMW1291",
                          "PMW1298",
                          "PMW1001"))

genotypes<-c("[ChrI-lacO, GBPD1-lacI::mCherry, SDC-1::GFP]",
             "[ChrI-lacO, lacI::mCherry, SDC-1::GFP]",
             "[ChrV-lacO, GBPD1-lacI::mCherry, SDC-1::GFP]")

p1<-ggplot(df[df$seqnames=="chrI",],aes(x=start/1e6,ymin=0,ymax=log2FoldChange)) +
  geom_linerange(color=ifelse(df$log2FoldChange[df$seqnames=="chrI"]>0,"royalblue","darkorange")) +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=targetdf$start[targetdf$seqnames=="chrI"]/1e6,linetype="dashed",color="grey")+
  geom_point(data=targetdf[targetdf$seqnames=="chrI",],aes(x=start/1e6,y=log2FoldChange),
             shape=15, color=c("firebrick3","firebrick3","black"),alpha=c(1,1,0), size=2) +
  geom_point(data=targetdf[targetdf$seqnames=="chrI",],aes(x=start/1e6,y=log2FoldChange+0.1),
             shape=15, color=c("forestgreen","black","black"),alpha=c(1,0,0), size=2) +
  facet_grid(group~.) +
  xlab("ChrI genomic position (Mb)") +
  ylab("Log<sub>2</sub>FC ( Bound / IPTG )") +
  coord_cartesian(ylim=ylim) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
  #geom_text(data=data.frame(x=7.38,y=0.5,label=genotypes,group=levels(df$group)),aes(x,y,label=label),
  #          inherit.aes=F,size=3,hjust=0)
p1

p2<-ggplot(df[df$seqnames=="chrV",],aes(x=start/1e6,ymin=0,ymax=log2FoldChange)) +
  geom_linerange(color=ifelse(df$log2FoldChange[df$seqnames=="chrV"]>0,"royalblue","darkorange")) +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=targetdf$start[targetdf$seqnames=="chrV"]/1e6,linetype="dashed",color="grey")+
  facet_grid(group~.) +
  geom_point(data=targetdf[targetdf$seqnames=="chrV",],aes(x=start/1e6,y=log2FoldChange),
             shape=15, color=c("black","black","firebrick3"),alpha=c(0,0,1),size=2) +
  geom_point(data=targetdf[targetdf$seqnames=="chrV",],aes(x=start/1e6,y=log2FoldChange+0.1),
             shape=15, color=c("black","black","forestgreen"),alpha=c(0,0,1), size=2) +
  xlab("ChrV genomic position (Mb)") +
  ylab("Log<sub>2</sub>FC ( Bound / IPTG )") +
  coord_cartesian(ylim=ylim) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
  #geom_text(data=data.frame(x=20.6,y=0.5,label=genotypes,group=levels(df$group)),aes(x,y,label=label),
  #          inherit.aes=F,size=3,hjust=0)
p2
p<-ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend="bottom")
p
ggsave(paste0(workDir,"/forPaper/","Fig_zoomed_",regionSize/1e6,"Mb_lineplots.pdf"),p,width=18,height=13,unit="cm")



# grFix=GRanges("chrI",Ranges(min(df[df$seqnames=="chrI","start"]),max(df[df$seqnames=="chrI","end"])))
#
# params <- pgParams(
#   chrom = as.character(seqnames(grFix)), chromstart = start(grFix), chromend = end(grFix),
#   assembly = "ce11",
#   just = c("left", "top"),
#   default.units = "cm",
#   x=0,
#   width = 8,
#   fontsize = 10,
#   spaceHeight=0, spaceWidth=0
# )
#
# pageCreate(width = 8, height = 12, params=params, showGuides=FALSE)
# p<-plotGenes(params=params,
#              y=9.5, height=1.5,
#              fontsize=8)
# annoGenomeLabel(
#   plot = p,
#   params=params, scale="Mb",
#   y = 11, height = 2
# )
#
#
#
# grFix=GRanges("chrV",IRanges(min(df[df$seqnames=="chrV","start"]),max(df[df$seqnames=="chrV","end"])))
#
# params <- pgParams(
#   chrom = as.character(seqnames(grFix)), chromstart = start(grFix), chromend = end(grFix),
#   assembly = "ce11",
#   just = c("left", "top"),
#   default.units = "cm",
#   x=1.4,
#   width = 8,
#   fontsize = 10,
#   spaceHeight=0, spaceWidth=0
# )
#
# pageCreate(width = 12, height = 12, params=params, showGuides=FALSE)
# pp<-plotGenes(params=params,
#              y=9.5, height=1.5,
#              fontsize=8)
# annoGenomeLabel(
#   plot = pp,
#   params=params, scale="Mb",
#   y = 11, height = 2
# )
#
# plotGG(p2, x=0, y=0, width=10.3, height=9,params=params)
