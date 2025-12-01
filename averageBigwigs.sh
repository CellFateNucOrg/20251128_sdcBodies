#! /usr/bin/bash
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=2


source $CONDA_ACTIVATE deeptools

workDir=/mnt/external.data/MeisterLab/jsemple/20251118_sdcBodies
runName=res01_allSamples_minAbund5_minSamples3
samplesheet=$workDir/fileList_extended.csv
sizeFactorFile=$( ls $workDir/$runName/other/deseq2/*.deseq2.sizefactors.tsv | head -1 )
blackListFile=$workDir/$runName/custom/tracks/c_elegans.PRJNA13758.WS295.blacklist.bed
outDir=$workDir/$runName/custom/tracks/bigwigs_fr

contrastsFile=$workDir/contrasts.csv
outDir2=$workDir/$runName/custom/tracks/bigwigs_compare

mkdir -p $outDir
mkdir -p $outDir2

calculate_inverses() {
    local numbers=("$@")
    local inverses=()

    for num in "${numbers[@]}"; do
        if [ "$num" != 0 ]; then
            inverse=$(echo "scale=10; 1 / $num" | bc -l)
            inverses+=("$inverse")
        else
            inverses+=("undefined")  # or you can use "" or 0 if that suits your case
        fi
    done
    echo "${inverses[@]}"
}

#### bigwigAverage ######

groups=(`tail -n +2 "$samplesheet" | cut -f1,6 -d"," | awk -F',' '{print $2}' | sort | uniq`)

for group in "${groups[@]}"; do
  echo "samples in group ${group}:"
  samples=(`awk -F',' -v val="$group" 'NR>1 && $6 == val { print $1 }' "$samplesheet"`)
  echo "${samples[@]}"

  sizeFactors=()
  for sample in "${samples[@]}"; do
    sizeFactors+=(`awk -F'\t' -v key="$sample" '$1 == key { print $2 }' "$sizeFactorFile"`)
    echo "$sample  with size factor ${sizeFactors[-1]}"
  done
  scaleFactors=($(calculate_inverses ${sizeFactors[@]}))

  bigwigs_f=( "${samples[@]/%/.forward.bigWig}" )
  bigwigs_r=( "${samples[@]/%/.reverse.bigWig}" )
  bigwigs_f=( "${bigwigs_f[@]/#/${workDir}/star_salmon/bigwig/}" )
  bigwigs_r=( "${bigwigs_r[@]/#/${workDir}/star_salmon/bigwig/}" )

  # NOTE: naming of star_salmon bigwigs with forward/reverse seems to be inverted. not sure if this is always true.
  # check on IGV and rename if necessary
  if [ ! -f "${outDir}/${group}.reverse.avr.bw" ]; then
    echo "averaging forward strand bigwigs"
    bigwigAverage -b ${bigwigs_f[@]} -o ${outDir}/${group}.reverse.avr.bw --scaleFactors $(IFS=:; echo "${scaleFactors[*]}") -bl $blackListFile
  fi

  if [ ! -f "${outDir}/${group}.forward.avr.bw" ]; then
    echo "averaging reverse strand bigwigs"
    bigwigAverage -b ${bigwigs_r[@]} -o ${outDir}/${group}.forward.avr.bw --scaleFactors $(IFS=:; echo "${scaleFactors[*]}") -bl $blackListFile
  fi

  if [ ! -f "${outDir}/${group}.average_fr.bw" ]; then
    echo "averaging forward and reverse strand bigwigs"
    bigwigAverage -b ${outDir}/${group}.forward.avr.bw ${outDir}/${group}.reverse.avr.bw -o ${outDir}/${group}.average_fr.bw
  fi
done

#### bigwigCompare ######

contrasts=(`tail -n +2 "$contrastsFile" | cut -f1 -d","`)
references=(`tail -n +2 "$contrastsFile" | cut -f3 -d","`)
targets=(`tail -n +2 "$contrastsFile" | cut -f4 -d","`)

for ((i=0; i<${#contrasts[@]}; i++)); do
    echo ${contrasts[i]} ${references[i]} ${targets[i]}
    if [ ! -f "${outDir2}/${contrasts[i]}.average_fr.log2fc.bw" ]; then
      bigwigCompare -b1 $outDir/${targets[i]}.average_fr.bw -b2 ${outDir}/${references[i]}.average_fr.bw -o ${outDir2}/${contrasts[i]}.average_fr.log2fc.bw
    fi
done

