#!/bin/bash

########################### BACKGROUND ###########################
# Intron clustering performed using the protocol specified in    #
# leafcutter `differntial splicing` instructions:                #
# https://davidaknowles.github.io/leafcutter/articles/Usage.html #
##################################################################

STAR_DIR=$1 #directory containing STAR alignments; for our purposes here, each bam possessed the following suffix: `*Aligned.sortedByCoord.out.bam`

#Estimate intron usage per sample (BAM) using regtools
for BAM in STAR_DIR
do
	#Variable prep
	library=`echo "${stardir}" | awk -F"/" '{print $(NF)}'`
	bamfile=$stardir/${library}Aligned.sortedByCoord.out.bam

	#Copy library to pwd
	echo "Copying $library"
	cp $bamfile .

	#Index library
	echo "Indexing $library"
	samtools index -@ 16 ${library}Aligned.sortedByCoord.out.bam

	#Generate regtools to quantify intron usage; parameters default as specified by leafcutter v.0.2.9
	echo "Generating ${library}.junc"
	regtools junctions extract -s 0 -a 8 -m 50 -M 500000 ${library}Aligned.sortedByCoord.out.bam -o ${library}Aligned.sortedByCoord.out.bam.junc

	#Remove intermediate files from pwd
	rm ${library}Aligned.sortedByCoord.out.bam
	rm ${library}Aligned.sortedByCoord.out.bam.bai

	#Write junc to txt file for clustering
	echo ${library}Aligned.sortedByCoord.out.bam.junc >> junc_files.txt
done

#Cluster intron usage estimates using leafcutter (again with default parameters specified in leafcutter v.0.2.9)
#Path to leafcutter executables can be modified below
python ~/bin/leafcutter/clustering/leafcutter_cluster_regtools.py -j junc_files.txt -m 50 -o All_samples -l 500000
