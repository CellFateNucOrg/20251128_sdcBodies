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

# which paramsheet goes with which contrasts file (relative to workDir); the contrasts file
# can be the csv or the pipeline YAML (e.g. with formula contrasts)
runSets <- list(
  list(params = "scripts/params.yml",        contrasts = "contrasts.csv"),
  list(params = "scripts/params_raptor.yml", contrasts = "contrasts_raptor.yml")
)


### functions ------
# contrast ids from a contrasts csv/tsv or YAML file, in file order
contrast_ids <- function(file) {
  if (grepl("\\.ya?ml$", file, ignore.case = TRUE)) {
    sapply(read_yaml(file)$contrasts, `[[`, "id")
  } else {
    read.csv(file, header = TRUE, sep = if (grepl("\\.tsv$", file)) "\t" else ",")$id
  }
}

annotate_results <- function(gtf, pattern="\\.deseq2\\.results\\.tsv$", path=".", runName=""){
  fileList<-list.files(paste0(path,"/tables/differential/",runName),pattern=pattern, full.names=T)
  if (length(fileList) == 0) warning(runName, ": no results files in ", paste0(path,"/tables/differential/",runName), call. = FALSE)
  outDir <- paste0(path,"/custom/txt/",runName)
  # remove annotated tables from earlier runs so stale files are not combined
  unlink(list.files(outDir, pattern="\\.deseq2\\.results_annotated\\.tsv$", full.names=T))
  for(f in fileList){
    print(f)
    df <- read.delim(f, header=T)
    if(!("gene_id" %in% colnames(df))){
      df$gene_id<-row.names(df)
      df$gene_name<-row.names(df)
    }
    # contrast id = file name without the "_study" (if present) and the results suffix
    file_name <- sub("(_study)?\\.deseq2\\.results\\.tsv$", "", basename(f))
    df <- inner_join(df, data.frame(gtf), by=c("gene_id" = "gene_id"))
    df$group <- file_name
    write.table(df, paste0(outDir,"/",file_name, ".deseq2.results_annotated.tsv"),
                row.names=F, quote=F, sep="\t", col.names=T)
  }
}


combine_results <- function(pattern="\\.deseq2\\.results_annotated\\.tsv$",path=".",runName=""){
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

# run names and the contrasts file of each run
runInfo <- do.call(rbind, lapply(runSets, function(s) data.frame(
  runName   = sapply(read_yaml(file.path(workDir, s$params)), `[[`, "paramset_name"),
  contrasts = file.path(workDir, s$contrasts))))
runInfo <- runInfo[!duplicated(runInfo$runName), ]


for(k in seq_len(nrow(runInfo))){
  runName <- runInfo$runName[k]
  message("== ", runName)
  ids <- contrast_ids(runInfo$contrasts[k])

  dir.create(paste0(workDir,"/diffabund/custom/rds/",runName), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(workDir,"/diffabund/custom/txt/",runName), showWarnings = FALSE, recursive = TRUE)

  ## annotate results tables -----
  annotate_results(gtf,path=paste0(workDir,"/diffabund"),runName=runName)

  results<-combine_results(path=paste0(workDir,"/diffabund"),runName=runName)
  if (is.null(results) || nrow(results) == 0) {
    warning(runName, ": no annotated results, nothing saved", call. = FALSE)
    next
  }
  results$seqnames<- factor(results$seqnames, levels=c("I", "II", "III", "IV", "V", "X", "MtDNA"))

  # contrasts in contrasts-file order; groups not in the file are kept (at the end) with a warning
  found <- unique(results$group)
  extra <- setdiff(found, ids)
  if (length(extra) > 0) warning(runName, ": groups not in ", basename(runInfo$contrasts[k]), ": ",
                                 paste(extra, collapse = ", "), call. = FALSE)
  results$group<-factor(results$group, levels=c(intersect(ids, found), extra))
  message("   groups: ", paste(levels(results$group), collapse = ", "))

  saveRDS(results, paste0(workDir,"/diffabund/custom/rds/",runName,"/results_annotated.RDS"))
}
