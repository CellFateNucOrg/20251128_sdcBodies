library(dplyr)
library(stringr)
library(yaml)

serverPath="/Volumes/meister.data"
workDir=paste0(serverPath,"/kalyan/sequencing/20261028_kg_ill_rna_l3_wPM1486_autosomal_smc1")
setwd(workDir)
fq<-read.csv(paste0(workDir,"/fileList.txt"),header=F)

nrow(fq)
idx<-seq(1, nrow(fq),by=2)

df<-data.frame(sample=NA,fastq_1=fq[idx,],fastq_2=fq[idx+1,],strandedness="auto")

df$sample<-sapply(strsplit(df$fastq_1,'/'), "[[", 8)
df$lane<-str_split_i(df$sample,"_",6)

df$sample<-paste0(gsub("_MKRN.*","",df$sample),"_",df$lane)

df$strain<-"PMW1486"
df$condition<-str_split_i(df$sample,"_",1)
df$replicate<-paste0("rep",gsub("S","",str_split_i(df$sample,"_",3)))
df$group<-paste0(df$condition,"_",df$strain)
df$lacO<-"oxTi481"

write.csv(df,file=paste0(workDir,"/samplesheet.csv"),quote=F,row.names=F)



contrasts<-data.frame(id=c("PWM1486_Bound_vs_IPTG"),
                      variable="group",
                      reference=c("IPTG_PMW1486"),
                      target=c("NO_PMW1486"),
                      blocking="replicate")

write.csv(contrasts,file=paste0(workDir,"/contrasts.csv"),quote=F,row.names=F)

q <- function(x) { attr(x, "quoted") <- TRUE; x }
contrastsyml<-list(contrasts=list(
  list(
    id = "PWM1486_Bound_vs_IPTG",
    formula = q("~ ageZnorm + replicate + group"),
    make_contrasts_str = q("groupNO_PMW1486")
    )
  )
)

write_yaml(contrastsyml, paste0(workDir,"/contrasts_raptor.yml"))
