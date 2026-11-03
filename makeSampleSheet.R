library(dplyr)
library(stringr)

workDir="/Volumes/meister.data/sinem/sequencing/20261028_se_ill_rna_l3_wPM1556_autosomal_bmm1"
fq<-read.csv(paste0(workDir,"/fileList.txt"),header=F)

nrow(fq)
idx<-seq(1, nrow(fq),by=2)

df<-data.frame(sample=NA,fastq_1=fq[idx,],fastq_2=fq[idx+1,],strandedness="unstranded")

df$sample<-sapply(strsplit(df$fastq_1,'/'), "[[", 8)
df$lane<-str_split_i(df$sample,"_",5)

df$sample<-paste0(gsub("_MKRN.*","",df$sample),"_",df$lane)

df$strain<-"PMW1556"
df$condition<-str_split_i(df$sample,"_",2)
df$replicate<-paste0("rep",gsub("S","",str_split_i(df$sample,"_",1)))
df$group<-paste0(df$condition,"_",df$strain)
df$lacO<-"oxTi481"

write.csv(df,file=paste0(workDir,"/samplesheet.csv"),quote=F,row.names=F)



contrasts<-data.frame(id=c("PWM1556_Bound_vs_IPTG"),
                      variable="group",
                      reference=c("IPTG_PMW1556"),
                      target=c("No_PMW1556"),
                      blocking="replicate")

write.csv(contrasts,file=paste0(workDir,"/contrasts.csv"),quote=F,row.names=F)

