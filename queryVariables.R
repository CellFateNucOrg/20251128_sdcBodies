


source("./variables.R")

print(paste(c(genomeVer,
        gsub("^\\/Volumes","/mnt",gtfFile),
        runName,
        gsub("^\\/Volumes","/mnt",sampleSheetFile),
        gsub("^\\/Volumes","/mnt",countsFile),
        gsub("^\\/Volumes","/mnt",lengthsFile),
        gsub("^\\/Volumes","/mnt",contrastsFile),
        ifelse(shrink,"true","false"),
        minAbund,
        minSamples),collapse=" "))
