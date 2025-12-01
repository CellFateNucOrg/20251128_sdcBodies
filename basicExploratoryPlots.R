# Plotting significant gene number and LFC by chromosome, chromosome region (arm/center),
# volcano plots and heatmaps

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(plotly)
library(ggrepel)
library(rstatix)
library(htmlwidgets)
library(RColorBrewer)
library(ComplexHeatmap)
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

serverPath="/Volumes/external.data/MeisterLab"
#serverPath="Z:/MeisterLab"
source("./functions_exploratory_analysis.R")

workDir=paste0(serverPath,"/jsemple/20251118_sdcBodies")
#runName="/res01_allSamples_minAbund5_minSamples3"
runName="/res02_no1298IPB2_minAbund5_minSamples3"
prefix=""

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)


setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/plots"), showWarnings = FALSE, recursive = TRUE)





## interactive Volcano plots -------
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
results$padj[is.na(results$padj)]<-1
dir.create(paste0(workDir,runName,"/custom/plots/volcano"), showWarnings = FALSE, recursive = TRUE)
#rockman<-readRDS(paste0(serverPath,"/publicData/Various/chrRegions_Rockman2009.RDS"))
gr<-tableToGranges(results,sort=FALSE)
seqlevelsStyle(gr)<-"UCSC"
#ol<-findOverlaps(resize(gr,width=1,fix="start"),rockman,ignore.strand=T)
#gr$chrRegionType[queryHits(ol)] <- paste0(rockman$chrRegionType[subjectHits(ol)])
#gr$chrRegion[queryHits(ol)] <- paste0(seqnames(rockman)[subjectHits(ol)],"_",rockman$chrRegionType[subjectHits(ol)])
res<-as.data.frame(gr)
contrasts

lfcVal=0.5
padjVal=0.05

# free scale
for(i in 1:length(contrasts$id)){
  tmp<-res[res$group==contrasts$id[i],c("padj","log2FoldChange","gene_name","gene_biotype","seqnames")]
  tmp$mlog10padj<--log10(tmp$padj)
  tmp$significant<-"NS"
  tmp$significant[tmp$log2FoldChange>lfcVal & tmp$padj<padjVal]<- "Up"
  tmp$significant[tmp$log2FoldChange< -lfcVal & tmp$padj<padjVal]<- "Down"

  tmp$significant<-factor(tmp$significant,levels=c("Up","Down","NS"))


  volcano_plot <- ggplot(tmp, aes(x = log2FoldChange, y = -log10(padj),
                                  color=significant,
                                  text = paste("Gene:", gene_name,
                                               "<br>Biotype:", gene_biotype,
                                               "<br>Chromosome:", seqnames,
                                               "<br>log2FC:", round(log2FoldChange, 2),
                                               "<br>padj:", signif(padj, 3)))) +
    geom_point(size=1.5,alpha=0.4) +
    scale_color_manual(values = c("Up" = "#E41A1C",
                                  "Down" = "#377EB8",
                                  "NS" = "lightgrey"),
                       guide = "none") +
    theme_minimal() + ggtitle(contrasts$id[i]) + xlab("Log2 Fold Change")

  volcano_plot1<-volcano_plot + geom_text_repel(
    data = subset(tmp, padj < 0.05 & abs(log2FoldChange) > 0.5),
    aes(label = gene_name),
    box.padding = 0,
    #max.overlaps = Inf,
    size = 2,
    force = 0.5, force_pull=1
  )
  volcano_plot1
  ggsave(filename = paste0(workDir,runName,"/custom/plots/volcano/volcano_",contrasts$id[i],"_free.png"),
         plot = volcano_plot1, width = 8, height = 6, dpi = 300,bg="white")
  volcano_plot
  interactive_volcano <- ggplotly(volcano_plot,tooltip="text")
  saveWidget(interactive_volcano, file = paste0(workDir,runName,"/custom/plots/volcano/volcano_",contrasts$id[i],"_free.html"))
}

# zoom LFC -5 to 5
for(i in 1:length(contrasts$id)){
  tmp<-res[res$group==contrasts$id[i],c("padj","log2FoldChange","gene_name","gene_biotype","seqnames")]
  tmp$mlog10padj<--log10(tmp$padj)
  tmp$significant<-"NS"
  tmp$significant[tmp$log2FoldChange>lfcVal & tmp$padj<padjVal]<- "Up"
  tmp$significant[tmp$log2FoldChange< -lfcVal & tmp$padj<padjVal]<- "Down"

  tmp$significant<-factor(tmp$significant,levels=c("Up","Down","NS"))


  volcano_plot <- ggplot(tmp, aes(x = log2FoldChange, y = -log10(padj),
                                  color=significant,
                                  text = paste("Gene:", gene_name,
                                               "<br>Biotype:", gene_biotype,
                                               "<br>Chromosome:", seqnames,
                                               "<br>log2FC:", round(log2FoldChange, 2),
                                               "<br>padj:", signif(padj, 3)))) +
    geom_point(size=1.5,alpha=0.4) +
    scale_color_manual(values = c("Up" = "#E41A1C",
                                  "Down" = "#377EB8",
                                  "NS" = "lightgrey"),
                       guide = "none") +
    theme_minimal() + ggtitle(contrasts$id[i]) + xlab("Log2 Fold Change") +
    coord_cartesian(xlim=c(-5,5))

  volcano_plot1<-volcano_plot + geom_text_repel(
    data = subset(tmp, padj < 0.05 & abs(log2FoldChange) > 0.5),
    aes(label = gene_name),
    box.padding = 0,
    #max.overlaps = Inf,
    size = 2,
    force = 0.5, force_pull=1
  )
  volcano_plot1
  ggsave(filename = paste0(workDir,runName,"/custom/plots/volcano/volcano_",contrasts$id[i],"_-5to5.png"),
         plot = volcano_plot1, width = 8, height = 6, dpi = 300,bg="white")
  volcano_plot
  interactive_volcano <- ggplotly(volcano_plot,tooltip="text")
  saveWidget(interactive_volcano, file = paste0(workDir, runName,"/custom/plots/volcano/volcano_",contrasts$id[i],"_-5to5.html"))
}


## Up down by chromosome  -----

dir.create(paste0(workDir, runName, "/custom/plots/byChromosome"), showWarnings = FALSE, recursive = TRUE)

lfcVal=0
padjVal=0.05

nuclear<-seqlevels(Celegans)[1:6]
resA<-res[res$seqnames %in% nuclear,]
resA$padj[is.na(resA$padj)]<-1
resA<-resA[resA$padj<padjVal & abs(resA$log2FoldChange)>lfcVal,]


### lfc -----
obsCounts<-data.frame(resA) %>% group_by(group,seqnames) %>%
  summarize(count = n())
wilcoxt<-data.frame(resA) %>% group_by(group) %>% wilcox_test(log2FoldChange~seqnames,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)
ylimits<-c(-3,3)
p<-ggplot(resA,aes(x=seqnames,y=log2FoldChange)) +
  geom_boxplot(aes(fill=seqnames),outlier.shape=NA) +
  facet_wrap(~group, scales="free") +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40") +
  scale_fill_brewer(palette="Blues")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.signif",y.position=ylimits[2]*0.95,
                     remove.bracket = T,angle=0,hide.ns=T) +
  ggtitle(paste0("Log2FoldChange for significant genes (padj<",padjVal," |LFC|>",lfcVal,")"))+
  theme(legend.position = "none")
p
ggsave(filename = paste0(workDir, runName, "/custom/plots/byChromosome/lfcByChr_",
                         contrasts$id[i], "_padj",padjVal,"_lfc",lfcVal,".png"),
       plot = p, width = 11, height = 8, dpi = 300,bg="white")


### counts -----
resA$upVdown<-NA
resA$upVdown[resA$log2FoldChange>lfcVal & resA$padj<padjVal]<-"up"
resA$upVdown[resA$log2FoldChange<lfcVal & resA$padj<padjVal]<-"down"
resA$upVdown<-factor(resA$upVdown, levels=c("up","down"))

obsCounts<-data.frame(resA) %>% group_by(group,seqnames,upVdown) %>%
  summarize(count = n())
obsCounts$y<-ifelse(obsCounts$upVdown=="up",0.95,0.05)
p<-ggplot(resA,aes(x=seqnames,fill=upVdown)) +
  geom_bar(position="fill") +
  facet_wrap(~group) +
  geom_hline(yintercept=0.5,linetype="dashed",color="black") +
  geom_text(data=obsCounts,aes(label=count,y=y),color="black",angle=0, size=3) +
  ggtitle(paste0("Fraction of up/down significant genes (padj<",padjVal," |LFC|>",lfcVal,")"))+
  scale_fill_brewer(palette = "Accent")
p
ggsave(filename = paste0(workDir, runName, "/custom/plots/byChromosome/countsByChr_",
                         contrasts$id[i],"_padj",padjVal,"_lfc",lfcVal,".png"),
       plot = p, width = 11, height = 8, dpi = 300,bg="white")




## Heatmap of significant genes -----
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir, runName, "/custom/plots/heatmaps"), showWarnings = FALSE, recursive = TRUE)


padjVal=0.05

makeHeatmapPlot(results, numSamplesSignificant=1, lfcVal=0, padjVal, direction="both",
                groupsToPlot=contrasts$id, setName="",
                chromosomes="all")


