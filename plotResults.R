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
source("./functions_exploratory_analysis.R")


contrasts<-read.csv(contrastsFile,sep=",",header=T)


results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))

shrink=ifelse(grepl("noShrink",runName),FALSE,TRUE)
if(shrink==TRUE){
  ylimits=c(-0.5,0.5)
} else {
  ylimits=c(-1.5,1.5)
}

dir.create(paste0(workDir, runName, "/custom/plots/byChromosome"), showWarnings = FALSE, recursive = TRUE)
df<-results[results$log2FoldChange!=0 & !is.na(results$padj),]
obsCounts<-data.frame(df) %>% group_by(group,seqnames) %>%
  summarize(count = n())
wilcoxt<-data.frame(df) %>% group_by(group) %>% wilcox_test(log2FoldChange~seqnames,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)# %>%

p<-ggplot(df, aes(x=seqnames, y=log2FoldChange)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~group) +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  labs(title=paste0("Differential expression by chromosome (",prefix,")"),
       x="Chromosome",
       y="Log2 Fold Change")  +
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]),color="blue",angle=90) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=ylimits[2]*0.20,
                     remove.bracket = T,angle=90,hide.ns=T,vjust=0.5)
p
ggsave(paste0(workDir,runName,"/custom/plots/byChromosome/",prefix,"_boxplot_byChr.pdf"), p, width=19, height=13, unit="cm")

# num significant
results %>% filter(!is.na(padj), padj < 0.05) %>% group_by(group,seqnames) %>%
  summarize(n_significant = n())


### chromosome_region ------
dir.create(paste0(workDir, runName, "/custom/plots/byChrRegion"), showWarnings = FALSE, recursive = TRUE)
gr<-tableToGranges(results)
seqlevelsStyle(gr) <- "UCSC"
chrRegions<-readRDS(paste0(serverPath,"/publicData/Various/Rockman_2009_armsVcenter/chrRegions_Rockman2009.RDS"))
chrRegions<-tableToGranges(chrRegions)
ol<-findOverlaps(gr,chrRegions,ignore.strand=T)
gr$chr_region <- NA
gr$chr_region[queryHits(ol)] <- chrRegions$name[subjectHits(ol)]

df<-as.data.frame(gr[gr$log2FoldChange!=0 & !is.na(gr$padj)])
obsCounts<-data.frame(df) %>% group_by(group,chr_region) %>%
  summarize(count = n())
wilcoxt<-data.frame(df) %>% group_by(group) %>% wilcox_test(log2FoldChange~chr_region,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

p<-ggplot(df, aes(x=chr_region, y=log2FoldChange)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~group) +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  labs(title=paste0("Differential expression by chromosome region type (",prefix,")"),
       x="Chromosome region type",
       y="Log2 Fold Change") +
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]),color="blue",angle=90) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=ylimits[2]*0.20,
                     remove.bracket = T,angle=90,hide.ns=T,vjust=0.5)
p
ggsave(paste0(workDir,runName,"/custom/plots/byChrRegion/",prefix,"_boxplot_byChrRegionType.pdf"), p, width=19, height=13, unit="cm")

data.frame(gr) %>% filter(!is.na(padj), padj < 0.2) %>% group_by(group,chr_region) %>%
  summarize(n_significant = n())


### anchors -------------
anchors<-import(paste0(serverPath,"/publicData/Various/Das_2024_anchors/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed"))
gr<-tableToGranges(results)
seqlevelsStyle(gr) <- "UCSC"
gr$anchor<-"autosome"
gr[seqnames(gr)=="chrX"]$anchor<-"nonAnchor"
ol<-findOverlaps(gr, anchors, ignore.strand=T)
gr$anchor[queryHits(ol)] <- "Anchor"
gr$anchor <-factor(gr$anchor,levels=c("autosome","nonAnchor","Anchor"))


df<-as.data.frame(gr[gr$log2FoldChange!=0 & !is.na(gr$padj)])
obsCounts<-data.frame(df) %>% group_by(group,anchor) %>%
  summarize(count = n())
wilcoxt<-data.frame(df) %>% group_by(group) %>% wilcox_test(log2FoldChange~anchor,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

p<-ggplot(df, aes(x=anchor, y=log2FoldChange)) +
  geom_boxplot(outlier.shape=NA) +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  labs(title=paste0("Differential expression, anchor vs non-anchor genes (",prefix,")"),
       x="Anchor Status",
       y="Log2 Fold Change") +
  facet_wrap(~group) +
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]),color="blue",angle=90) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=ylimits[2]*0.20,
                     remove.bracket = T,angle=90,hide.ns=T,hjust=0.5)
p
ggsave(paste0(workDir,runName,"/custom/plots/byChrRegion/",prefix,"_boxplot_anchorVnotanchor.pdf"), p, width=19, height=13, unit="cm")


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

# chrI
minBin=min(as.numeric(df$binNum[df$seqnames=="chrI"]))
maxBin=max(as.numeric(df$binNum[df$seqnames=="chrI"]))

subdf<-df[df$seqnames=="chrI" & df$group=="PMW1001_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

if(shrink){
  ylim=ylimits
} else {
  ylim=ylimits
}


p1<-ggplot(subdf, aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target),outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("ChrI ",binWidth/1e3,"kb bins, chrV_GBPD1-lacI_BoundvIPTG")) +
  geom_hline(yintercept=0,alpha=0.3) +
  scale_x_discrete(name="ChrI (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(minBin,maxBin,5)*binWidth/1e6) +
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90, size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)
p1

subdf<-df[df$seqnames=="chrI" & df$group=="PMW1291_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

if(shrink){
  ylim=ylimits
} else {
  ylim=ylimits
}
p1a<-ggplot(subdf,aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target),outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("ChrI ",binWidth/1e3,"kb bins, chrI_GBPD1-lacI_BoundvIPTG")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrI (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(minBin,maxBin,5)*binWidth/1e6) +
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90,size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)
p1a

subdf<-df[df$seqnames=="chrI" & df$group=="PMW1298_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

if(shrink){
  ylim=ylimits
} else {
  ylim=ylimits
}

p1b<-ggplot(subdf,aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target),outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("ChrI ",binWidth/1e3,"kb bins, chrI_lacI_BoundvIPTG")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrI (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(minBin,maxBin,5)*binWidth/1e6) +
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90, size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)
p1b

p<-ggarrange(p1,p1a,p1b,nrow=3)
p

ggsave(paste0(workDir,runName,"/custom/plots/byChrRegion/",prefix,"_chrI_",binWidth/1e3,"kb.pdf"),p,width=32,height=21,unit="cm")


# chrV
if(shrink==TRUE){
  ylimits=c(-0.05,0.05)
} else {
  ylimits=c(-1.5,1.5)
}

minBin=min(as.numeric(df$binNum[df$seqnames=="chrV"]))
maxBin=max(as.numeric(df$binNum[df$seqnames=="chrV"]))

subdf<-df[df$seqnames=="chrV" & df$group=="PMW1001_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

if(shrink){
  #ylim=c(ylimits[1]*6/5,ylimits[2]*2/5)
  ylim=ylimits
} else {
  ylim=ylimits
}
p2<-ggplot(subdf, aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target),outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("ChrV ",binWidth/1e3,"kb bins, chrV_GBPD1-lacI_BoundvIPTG")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrV (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(1,maxBin-minBin-1,5)*binWidth/1e6)+
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90,size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)
p2

subdf<-df[df$seqnames=="chrV" & df$group=="PMW1291_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

if(shrink){
  ylim=ylimits
} else {
  ylim=ylimits
}

p2a<-ggplot(subdf, aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target), outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("ChrV ",binWidth/1e3,"kb bins, chrI_GBPD1-lacI_BoundvIPTG")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrV (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(1,maxBin-minBin-1,5)*binWidth/1e6)+
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90, size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)
p2a


subdf<-df[df$seqnames=="chrV" & df$group=="PMW1298_Bound_vs_IPTG",]
subdf<-subdf[subdf$log2FoldChange!=0 & !is.na(subdf$padj),]
obsCounts<-data.frame(subdf) %>% group_by(group,binNum) %>%
  summarize(count = n())
wilcoxt<-data.frame(subdf) %>% group_by(group) %>% wilcox_test(log2FoldChange~binNum,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)

if(shrink){
  ylim=ylimits
} else {
  ylim=ylimits
}

p2b<-ggplot(subdf,aes(x=binNum,y=log2FoldChange)) +
  geom_boxplot(aes(fill=target), outliers=F,show.legend=F) +
  coord_cartesian(ylim=ylim) +
  scale_fill_manual(values = c("grey90", "red"))+
  ggtitle(paste0("ChrV ",binWidth/1e3,"kb bins, chrI_lacI_BoundvIPTG")) +
  geom_hline(yintercept=0,alpha=0.3)+
  scale_x_discrete(name="ChrV (Mb)",breaks=seq(minBin,maxBin,5),labels=seq(1,maxBin-minBin-1,5)*binWidth/1e6)+
  geom_text(data=obsCounts,aes(label=count,y=ylim[1]*0.9),color="blue",angle=90, size=4) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.format",y.position=(ylim[2])*0.5,
                     remove.bracket = T,angle=90,hide.ns=T,size=3)
p2b

p<-ggarrange(p2,p2a,p2b,nrow=3)
p
ggsave(paste0(workDir,runName,"/custom/plots/byChrRegion/",prefix,"_chrV_",binWidth/1e3,"kb.pdf"),p,width=32,height=21,unit="cm")



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


p1<-ggplot(df[df$seqnames=="chrI",],aes(x=start/1e6,ymin=0,ymax=log2FoldChange)) +
  geom_linerange(color=ifelse(df$log2FoldChange[df$seqnames=="chrI"]>0,"royalblue","darkorange")) +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=targetdf$start[targetdf$seqnames=="chrI"]/1e6,linetype="dashed",color="grey")+
  theme_classic() +
  geom_point(data=targetdf[targetdf$seqnames=="chrI",],aes(x=start/1e6,y=log2FoldChange),
             shape=15, color=c("black","firebrick3","firebrick3"),alpha=c(0,1,1), size=2) +
  geom_point(data=targetdf[targetdf$seqnames=="chrI",],aes(x=start/1e6,y=log2FoldChange+0.07),
             shape=15, color=c("black","forestgreen","black"),alpha=c(0,1,0), size=2) +
  facet_grid(group~.) +
  xlab("ChrI genomic position (Mb)")
p1

p2<-ggplot(df[df$seqnames=="chrV",],aes(x=start/1e6,ymin=0,ymax=log2FoldChange)) +
  geom_linerange(color=ifelse(df$log2FoldChange[df$seqnames=="chrV"]>0,"royalblue","darkorange")) +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=targetdf$start[targetdf$seqnames=="chrV"]/1e6,linetype="dashed",color="grey")+
  theme_classic() +
  facet_grid(group~.) +
  geom_point(data=targetdf[targetdf$seqnames=="chrV",],aes(x=start/1e6,y=log2FoldChange),
             shape=15, color=c("firebrick3","black","black"),alpha=c(1,0,0),size=2) +
  geom_point(data=targetdf[targetdf$seqnames=="chrV",],aes(x=start/1e6,y=log2FoldChange+0.12),
             shape=15, color=c("forestgreen","black","black"),alpha=c(1,0,0), size=2) +
  xlab("ChrV genomic position (Mb)")
p2
p<-ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend="bottom")
ggsave(paste0(workDir,runName,"/custom/plots/byChrRegion/",prefix,"_zoomed_",regionSize/1e6,"Mb_lineplots.pdf"),p,width=19,height=19,unit="cm")
