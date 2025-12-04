library(dplyr)

source("./variables.R")

fq<-read.csv("./fileList.txt",header=F)

nrow(fq)
idx<-seq(1, nrow(fq),by=2)

df<-data.frame(sample=NA,fastq_1=fq[idx,],fastq_2=fq[idx+1,],strandedness="auto")

df$sample<-sapply(strsplit(df$fastq_1,'/'), "[[", 8)

df$sample<-gsub("W01","PMW1001",df$sample)
df$sample<-gsub("W91","PMW1291",df$sample)
df$sample<-gsub("W98","PMW1298",df$sample)
df$sample<-gsub("_IP_","_IPTG_",df$sample)
df$sample<-gsub("_NO_","_Bound_",df$sample)

df$sample

write.csv(df,file=paste0(workDir,"/fileList.csv"),quote=F,row.names=F)

samples="allSamples"
df$strain<-sapply(strsplit(as.character(df$sample),"_"), "[[", 1)
df$condition<-sapply(strsplit(as.character(df$sample),"_"), "[[", 2)
df$group<-paste0(df$strain,"_",df$condition)
df$replicate<-sapply(strsplit(as.character(df$sample),"_"), "[[", 3)
df
write.csv(df,file=paste0(workDir,"/fileList_",samples,".csv"),quote=F,row.names=F)


contrasts<-data.frame(id=c("PMW1001_Bound_vs_IPTG",
                             "PMW1291_Bound_vs_IPTG",
                             "PMW1298_Bound_vs_IPTG"),
                      variable="group",
                      reference=c("PMW1001_IPTG",
                                  "PMW1291_IPTG",
                                  "PMW1298_IPTG"),
                      target=c("PMW1001_Bound",
                               "PMW1291_Bound",
                               "PMW1298_Bound"),
                      blocking="replicate")

write.csv(contrasts,file=paste0(workDir,"/contrasts.csv"),quote=F,row.names=F)


# removing 1298_IP_B2 -----
samples="no1298IPB2"
ss<-read.csv(paste0(workDir,"/fileList_allSamples.csv"))
ss<-ss[ss$sample!="PMW1298_IPTG_B2",]
write.csv(ss,file=paste0(workDir,"/fileList_",samples,".csv"),quote=F,row.names=F)

## NOTE: as long as sample sheet has correct samples, you do not need to filter count/length/tpm matrix columns



# combine samples sheets for bulk and low input ------
## NOTE: need to add columns from other samples sheet and also add columns to counts/lengths/tpm matrices
samples="no1298IPB2_lowInput"
bulk<-read.csv(paste0(workDir,"/fileList_no1298IPB2.csv"))
bulk$seqType<-"bulk"
lowinput<-read.csv(paste0(serverPath,"/RNA_seq_BCN/202501_PM/Sinem/samplesheet_extended_compareBulk.csv"))

#lowinput$sample<-gsub("PWM","PMW",lowinput$sample)
#lowinput$strain<-gsub("PWM","PMW",lowinput$strain)
#lowinput$group<-gsub("PWM","PMW",lowinput$group)
#write.csv(lowinput,paste0(serverPath,"/RNA_seq_BCN/202501_PM/Sinem/samplesheet_extended_compareBulk.csv"),row.names=F,quote=F)

toRemove<-colnames(lowinput)[! (colnames(lowinput) %in% colnames(bulk))] # remove gene id columns for merging
lowinput[,toRemove]<-NULL
colOrder<-match(colnames(bulk),colnames(lowinput))

newss<-rbind(bulk, lowinput[,colOrder])
write.csv(newss,file=paste0("fileList_",samples,".csv"),quote=F,row.names=F)
print(paste0("writing ",workDir,"/fileList_",samples,".csv"))

# combine counts
bulk_cnts<-read.csv(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"),sep="\t")

lowinput_cnts<-read.csv(paste0(serverPath,"/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_counts.tsv"),sep="\t")

#colnames(lowinput_cnts)<-gsub("PWM","PMW",colnames(lowinput_cnts))
#write.table(lowinput_cnts,paste0(serverPath,"/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_counts.tsv"),sep="\t",row.names=F,quote=F)

keep<-lowinput_cnts$gene_id %in% bulk_cnts$gene_id
lowinput_cnts_filt<-lowinput_cnts[keep,]

cnts<-cbind(bulk_cnts,lowinput_cnts_filt[3:ncol(lowinput_cnts_filt)])
dim(cnts)
write.table(cnts,file=paste0(workDir,"/star_salmon/salmon.merged.gene_counts.",samples,".tsv"),sep="\t",row.names=F,quote=F)


# combine gene lengths
bulk_len<-read.csv(paste0(workDir,"/star_salmon/salmon.merged.gene_lengths.tsv"),sep="\t")

lowinput_len<-read.csv(paste0(serverPath,"/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_lengths.tsv"),sep="\t")

#colnames(lowinput_len)<-gsub("PWM","PMW",colnames(lowinput_len))
#write.table(lowinput_len,paste0(serverPath,"/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_lengths.tsv"),sep="\t",row.names=F,quote=F)


keep<-lowinput_len$gene_id %in% bulk_len$gene_id
lowinput_len_filt<-lowinput_len[keep,]

len<-cbind(bulk_len,lowinput_len_filt[3:ncol(lowinput_len_filt)])
dim(len)
write.table(len,file=paste0(workDir,"/star_salmon/salmon.merged.gene_lengths.",samples,".tsv"),sep="\t",row.names=F,quote=F)



# combine gene tpm (for raptor)
bulk_tpm<-read.csv(paste0(workDir,"/star_salmon/salmon.merged.gene_tpm.tsv"),sep="\t")

lowinput_tpm<-read.csv(paste0(serverPath,"/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_tpm.tsv"),sep="\t")

#colnames(lowinput_tpm)<-gsub("PWM","PMW",colnames(lowinput_tpm))
#write.table(lowinput_tpm,paste0(serverPath,"/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_tpm.tsv"),sep="\t",row.names=F,quote=F)

keep<-lowinput_tpm$gene_id %in% bulk_tpm$gene_id
lowinput_tpm_filt<-lowinput_tpm[keep,]

tpm<-cbind(bulk_tpm,lowinput_tpm_filt[3:ncol(lowinput_tpm_filt)])

write.table(tpm,file=paste0(workDir,"/star_salmon/salmon.merged.gene_tpm.",samples,".tsv"),sep="\t",row.names=F,quote=F)

