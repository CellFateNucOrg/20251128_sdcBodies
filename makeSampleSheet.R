library(dplyr)
library(stringr)
library(yaml)

serverPath="/Volumes/meister.data"
workDir=paste0(serverPath,"/sinem/sequencing/20251118_Sinem_sdcBodies")
setwd(workDir)
fq<-read.csv(paste0(workDir,"/fileList.txt"),header=F)

nrow(fq)
idx<-seq(1, nrow(fq),by=2)

df<-data.frame(sample=NA,fastq_1=fq[idx,],fastq_2=fq[idx+1,],strandedness="auto")

df$sample<-sapply(strsplit(df$fastq_1,'/'), "[[", 8)
#df$lane<-str_split_i(basename(df$fastq_1),"_",6)

#df$sample<-paste0(gsub("_MKRN.*","",df$sample),"_",df$lane)

df$strain<-str_split_i(df$sample,"_",1)
df$strain<-gsub("W01","PMW1001",df$strain)
df$strain<-gsub("W91","PMW1291",df$strain)
df$strain<-gsub("W98","PMW1298",df$strain)
df$condition<-str_split_i(df$sample,"_",2)
df$replicate<-paste0("rep",gsub("S","",str_split_i(df$sample,"_",3)))
df$group<-paste0(df$condition,"_",df$strain)
df$lacO<-"oxTi481"
df$lacO[df$strain=="PMW1001"]<-"ttTi9115"

write.csv(df,file=paste0(workDir,"/samplesheet.csv"),quote=F,row.names=F)



contrasts<-data.frame(id=c("PWM1001_Bound_vs_IPTG","PWM1291_Bound_vs_IPTG","PWM1298_Bound_vs_IPTG"),
                      variable="group",
                      reference=c("IP_PMW1001","IP_PMW1291","IP_PMW1298"),
                      target=c("NO_PMW1001","NO_PMW1291","NO_PMW1298"),
                      blocking="replicate")

write.csv(contrasts,file=paste0(workDir,"/contrasts.csv"),quote=F,row.names=F)

q <- function(x) { attr(x, "quoted") <- TRUE; x }
contrastsyml<-list(contrasts=list(
  list(
    id = "PWM1001_Bound_vs_IPTG",
    formula = q("~ ageZnorm + replicate + group"),
    make_contrasts_str = q("groupNO_PMW1001")
    ),
  list(
    id = "PWM1291_Bound_vs_IPTG",
    formula = q("~ ageZnorm + replicate + group"),
    make_contrasts_str = q("groupNO_PMW1291")
  ),
  list(
    id = "PWM1298_Bound_vs_IPTG",
    formula = q("~ ageZnorm + replicate + group"),
    make_contrasts_str = q("groupNO_PMW1298")
  )
  )
)

write_yaml(contrastsyml, paste0(workDir,"/contrasts_raptor.yml"))
