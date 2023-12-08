#!/usr/bin/bash

#==========#
# Get args #
#==========#

fastqListFile=$1 # Three-column tab-separated file with libraryID, path to pair1 reads, path to pair2 reads
gencodeTranscriptsFasta=$2 #gencode.v38.transcripts.fa.gz
outIndexDir=$3
outSalmonDir=$4
threads=8


#=====================#
# Create salmon index #
#=====================#

salmon index -t $gencodeTranscriptsFasta -i $outIndexDir --gencode


#=========================#
# Create output directory #
#=========================#

if [[ ! -d $outSalmonDir ]]; then
	mkdir -p $outSalmonDir
fi


#============#
# Run salmon #
#============#

while IFS=$'\t' read -r -a lineArray; do

	libraryID=${lineArray[0]}
	pair1_fastq=${lineArray[1]}
	pair2_fastq=${lineArray[2]}

	salmon quant \
	-i $outIndexDir \
	-l IU \
	-1 $pair1_fastq \
	-2 $pair2_fastq \
	-p $threads \
	-o $outSalmonDir/$libraryID

done < $fastqListFile
