suppressMessages(library("log4r"))

resultsRun="res03a"

samples="allSamples"
#samples="no1298IPB2"
#samples="no1298IPB2_lowInput"


#minAbund=5
minAbund=10

#minSamples=9
minSamples=16

filterNoncoding=T
filterMitochondrial=T

filterOscillating=T
raptor=F
if(raptor){
  filterOscillating=F
}

shrink=F

if(Sys.info()['sysname']=="Darwin"){
  serverPath="/Volumes/external.data/MeisterLab"
} else if (Sys.info()['sysname']=="Linux"){
  serverPath="/mnt/external.data/MeisterLab"
} else if (Sys.info()['sysname']=="Windows"){
  serverPath="Z:/MeisterLab"
}


workDir=paste0(serverPath,"/jsemple/20251118_sdcBodies")


runName=paste0("/",resultsRun,"_",samples,"_minAbund",minAbund,"_minSamples",
               minSamples, ifelse(filterNoncoding,"_pc",""),
               ifelse(filterOscillating,"_noOsc",""),
               ifelse(filterMitochondrial,"_noMito",""),
               ifelse(raptor,"_raptor",""),
               ifelse(shrink,"","_noShrink"))

prefix=""

dir.create(paste0(workDir,runName,"/logs"), showWarnings = FALSE, recursive = TRUE)

log_file<-paste0(workDir,runName,"/logs/variables.log")

file_logger <- logger("INFO", appenders = file_appender(log_file,append=F))
info(file_logger, paste0("workDir: ",workDir))
info(file_logger, paste0("runName: ",runName))
#file_appender(log_file, append = TRUE, layout = default_log_layout())


## sampleSheet
sampleSheetFile<-paste0(workDir,"/fileList_",samples,".csv")
info(file_logger, paste0("sampleSheet: ",sampleSheetFile))

## contrasts file
contrastsFile<-paste0(workDir,"/contrasts.csv")
info(file_logger, paste0("contrastsFile: ",contrastsFile))

## gtf file
genomeVer="WS295"
if(filterNoncoding & filterMitochondrial){
  gtfFile<-paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.protein_coding.noMito.gtf")
} else if (filterNoncoding){
  gtfFile<-paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.protein_coding.gtf")
} else {
  gtfFile<-paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf")
}
log4r::info(file_logger, paste0("gtf file: ",gtfFile))

## counts file
countsFile<-paste0(workDir,"/star_salmon/salmon.merged.gene_counts",
                               ifelse(filterNoncoding,".protein_coding",""),
                               ifelse(filterMitochondrial,".noMito",""),
                               ifelse(filterOscillating,".noOsc",""),
                               ifelse(samples=="allSamples","",paste0(".",samples)),
                               ".tsv")
log4r::info(file_logger, paste0("counts file: ",countsFile))

## lengths file
if(samples=="no1298IPB2_lowInput"){
  lengthsFile<-paste0(workDir,"/star_salmon/salmon.merged.gene_lengths.",samples,".tsv")
} else {
  lengthsFile<-paste0(workDir,"/star_salmon/salmon.merged.gene_lengths.tsv")
}
log4r::info(file_logger, paste0("lengths file: ",lengthsFile))


## tpm file
if(samples=="no1298IPB2_lowInput"){
  tpmFile<-paste0(workDir,"/star_salmon/salmon.merged.gene_tpm.", samples,".tsv")
} else {
  tpmFile<-paste0(workDir,"/star_salmon/salmon.merged.gene_tpm.tsv")
}

log4r::info(file_logger, paste0("tpm file: ",tpmFile))


if(Sys.info()['sysname']=="Darwin" |  Sys.info()['sysname']=="Windows"){
  devtools::session_info(to_file=paste0(workDir,runName,"/logs/variables_sessionInfo.txt"))
}

