library(dplyr)

fq<-read.csv("./fileList.txt",header=F)

nrow(fq)
idx<-seq(1, nrow(fq),by=2)

df<-data.frame(sample=NA,fastq_1=fq[idx,],fastq_2=fq[idx+1,],strandedness="auto")

df$sample<-sapply(strsplit(df$fastq_1,'/'), "[[", 8)

df$sample<-gsub("W01","PWM1001",df$sample)
df$sample<-gsub("W91","PWM1291",df$sample)
df$sample<-gsub("W98","PWM1298",df$sample)
df$sample<-gsub("_IP_","_IPTG_",df$sample)
df$sample<-gsub("_NO_","_Bound_",df$sample)

  df$sample

write.csv(df,file="fileList.csv",quote=F,row.names=F)


df$strain<-sapply(strsplit(as.character(df$sample),"_"), "[[", 1)
df$condition<-sapply(strsplit(as.character(df$sample),"_"), "[[", 2)
df$group<-paste0(df$strain,"_",df$condition)
df$replicate<-sapply(strsplit(as.character(df$sample),"_"), "[[", 3)
df

write.csv(df,file="fileList_extended.csv",quote=F,row.names=F)



contrasts<-data.frame(id=c("PWM1001_Bound_vs_IPTG",
                             "PWM1291_Bound_vs_IPTG",
                             "PWM1298_Bound_vs_IPTG"),
                      variable="group",
                      reference=c("PWM1001_IPTG",
                                  "PWM1291_IPTG",
                                  "PWM1298_IPTG"),
                      target=c("PWM1001_Bound",
                               "PWM1291_Bound",
                               "PWM1298_Bound"),
                      blocking="replicate")

write.csv(contrasts,file="contrasts.csv",quote=F,row.names=F)

# combine samples sheets
bulk<-read.csv(file="fileList_extended.csv")
bulk$seqType<-"bulk"
lowinput<-read.csv("/Volumes/external.data/MeisterLab/RNA_seq_BCN/202501_PM/Sinem/samplesheet_extended_compareBulk.csv")


toRemove<-colnames(lowinput)[! (colnames(lowinput) %in% colnames(bulk))]
lowinput[,toRemove]<-NULL
colOrder<-match(colnames(bulk),colnames(lowinput))

newss<-rbind(bulk, lowinput[,colOrder])
write.csv(newss,file="fileList_extended_compareLI.csv",quote=F,row.names=F)


# add blocking variable to contrasts
contrasts<-read.csv(file="contrasts.csv")
contrasts$blocking<-paste0("seqType;",contrasts$blocking)
write.csv(contrasts,file="contrasts_compareLI.csv",quote=F,row.names=F)

# combine counts
bulk_len<-read.csv("./star_salmon/salmon.merged.gene_counts.tsv",sep="\t")

lowinput_cnts<-read.csv("/Volumes/external.data/MeisterLab/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_counts.tsv",sep="\t")

keep<-lowinput_cnts$gene_id %in% bulk_cnts$gene_id
lowinput_cnts_filt<-lowinput_cnts[keep,]

cnts<-cbind(bulk_cnts,lowinput_cnts_filt[3:ncol(lowinput_cnts_filt)])

write.table(cnts,file="./star_salmon/salmon.merged.gene_counts.bulk_plus_lowinput.tsv",sep="\t",row.names=F,quote=F)



# combine gene lengths
bulk_len<-read.csv("./star_salmon/salmon.merged.gene_lengths.tsv",sep="\t")

lowinput_len<-read.csv("/Volumes/external.data/MeisterLab/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_lengths.tsv",sep="\t")

keep<-lowinput_len$gene_id %in% bulk_len$gene_id
lowinput_len_filt<-lowinput_len[keep,]

len<-cbind(bulk_len,lowinput_len_filt[3:ncol(lowinput_len_filt)])

write.table(len,file="./star_salmon/salmon.merged.gene_lengths.bulk_plus_lowinput.tsv",sep="\t",row.names=F,quote=F)



# combine gene tpm (for raptor)
bulk_tpm<-read.csv("./star_salmon/salmon.merged.gene_tpm.tsv",sep="\t")

lowinput_tpm<-read.csv("/Volumes/external.data/MeisterLab/RNA_seq_BCN/202501_PM/Sinem/star_salmon/salmon.merged.gene_tpm.tsv",sep="\t")

keep<-lowinput_tpm$gene_id %in% bulk_tpm$gene_id
lowinput_tpm_filt<-lowinput_tpm[keep,]

tpm<-cbind(bulk_tpm,lowinput_tpm_filt[3:ncol(lowinput_tpm_filt)])

write.table(tpm,file="./star_salmon/salmon.merged.gene_tpm.bulk_plus_lowinput.tsv",sep="\t",row.names=F,quote=F)
