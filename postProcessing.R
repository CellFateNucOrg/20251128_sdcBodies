## add annotation to results tables with gene names and locations.
## Combine all results tables into a single .rds object
## extract table of Number of significant up/down regulated genes by different thresholds

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(yaml)



serverPath="/Volumes/meister.data"
#serverPath="Z:/MeisterLab"

workDir=paste0(serverPath,"/sinem/sequencing/20261028_se_ill_rna_l3_wPM1556_autosomal_bmm1")

setwd(workDir)

params <- read_yaml(paste0(workDir,"/scripts/params.yml"))

runNames  <- sapply(params, `[[`, "paramset_name")

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)


### functions ------
annotate_results <- function(gtf, pattern="\\.deseq2\\.results\\.tsv", path=".", runName=""){
  fileList<-list.files(paste0(path,"/tables/differential/",runName),pattern=pattern, full.names=T)
  for(f in fileList){
    print(f)
    df <- read.delim(f, header=T)
    if(!("gene_id" %in% colnames(df))){
      df$gene_id<-row.names(df)
      df$gene_name<-row.names(df)
    }
    file_name <- gsub("_study\\.deseq2\\.results\\.tsv", "", basename(f))
    #group_name <- gsub(paste0(prefix,"_"),"",file_name)
    #print(group_name)
    df <- inner_join(df, data.frame(gtf), by=c("gene_id" = "gene_id"))
    df$group <- file_name #group_name
    write.table(df, paste0(path,"/custom/txt/",runName,"/",file_name, ".deseq2.results_annotated.tsv"),
                row.names=F, quote=F, sep="\t", col.names=T)
  }
}


combine_results <- function(pattern="\\.deseq2\\.results_annotated\\.tsv",path=".",runName=""){
  fileList <- list.files(paste0(path,"/custom/txt/",runName), pattern=pattern, full.names=T)
  results <- do.call(rbind, lapply(fileList, read.delim, header=T))
  return(results)
}

#-----------------
genomeVer<-"WS298"
gtf<-import(paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))
gtf <- gtf[gtf$type == "gene"]
mcols(gtf)<-mcols(gtf)[c("source","type","gene_id","gene_biotype","gene_name")]
gtf$source<-genomeVer
gtf<-sort(gtf)

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)


for(runName in runNames){
  dir.create(paste0(workDir,"/diffabund/custom/rds/",runName), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(workDir,"/diffabund/custom/txt/",runName), showWarnings = FALSE, recursive = TRUE)

  ## annotate results tables -----
  annotate_results(gtf,path=paste0(workDir,"/diffabund"),pattern="\\.deseq2\\.results\\.tsv",runName=runName)

  results<-combine_results(path=paste0(workDir,"/diffabund"),pattern="\\.deseq2\\.results_annotated\\.tsv",runName=runName)
  results$seqnames<- factor(results$seqnames, levels=c("I", "II", "III", "IV", "V", "X", "MtDNA"))

  results$group<-factor(results$group, levels=contrasts$id)

  # if (file.exists(paste0(workDir,"/diffabund/custom/rds/",runName,"/results_annotated.RDS"))) {
  #   response <- readline(prompt = paste(runName, "file already exists. Overwrite? [y/n]: "))
  #   if (tolower(response) == "y") {
  #     saveRDS(results, paste0(workDir,"/diffabund/custom/rds/",runName,"/results_annotated.RDS"))
  #   } else {
  #     print("file not overwritten")
  #   }
  # } else {
  saveRDS(results, paste0(workDir,"/diffabund/custom/rds/",runName,"/results_annotated.RDS"))
  #}
}


