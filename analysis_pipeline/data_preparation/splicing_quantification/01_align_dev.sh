#!/bin/bash

#==========#
# Get args #
#==========#

workDir=$1     # Working directory where outputs will be generated
metaData=$2    # The path to the `updated_metadata.tsv` file, provided in this directory
sampleID=$3    # The identifier of a sample to align, such as `GM18861` or `HG00148`
               # replicates are denoted based on batch and replicate number, such as `GM20815-13-1` and `GM20815-13-2`
threads=$4     # The number of threads to use for computation (set to 24 for our analysis)
kgpVCFdir=$5   # Directory containing per-chromosome 1000 Genomes VCF files with names such as 
               # `CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.eagle2-phased.v2.vcf.gz`
starIndex=$6   # Path to the STAR index directory (generated from `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`)
fastqDir=$7    # Path to directory containing FASTQ files, organized into subdirectories based on sequencing batch

#===================================#
# Set variables by parsing metadata #
#===================================#

coriell_id=`echo ${sampleID} | cut -f1 -d '-'`
participant_id=`cut -f1,3- ${metaData} | grep -P "${sampleID}\t" | cut -f2`
batch="batch"`cut -f1,3- ${metaData} | grep -P "${sampleID}\t" | cut -f3`
sex=`cut -f1,3- ${metaData} | grep -P "${sampleID}\t" | cut -f7`

cd ${workDir}
mkdir output
mkdir output/${sampleID}
mkdir vcf

## generate individual sample-specific VCF files with heterozygous variants
j=1
if [ ${sex} = "female" ]
then
    (i="X"
     echo "Processing chromosome ${i}."
     vcf_file=${kgpVCFdir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.eagle2-phased.v2.vcf.gz
     bcftools view --threads ${threads} --no-update -s ${participant_id} -v snps ${vcf_file} | bcftools view --threads ${threads} --no-update -e 'GT=".|."' -Oz -o vcf/${sampleID}.chr${i}.snps.vcf.gz
     bcftools view --threads ${threads} --no-update -i 'GT="het"' vcf/${sampleID}.chr${i}.snps.vcf.gz | bcftools norm -m+ | bcftools view --threads ${threads} -m2 -M2 -Oz -o vcf/${sampleID}.chr${i}.snps.het.vcf.gz
     tabix vcf/${sampleID}.chr${i}.snps.het.vcf.gz) &
    pids[${j}]=$!
    j=$((j+1))
fi
for i in {1..22}; do
    (echo "Processing chromosome ${i}."
     vcf_file=${kgpVCFdir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz
     bcftools view --threads 1 --no-update -s ${participant_id} -v snps ${vcf_file} | bcftools view --threads 1 --no-update -e 'GT=".|."' -Oz -o vcf/${sampleID}.chr${i}.snps.vcf.gz
     bcftools view --threads 1 --no-update -i 'GT="het"' vcf/${sampleID}.chr${i}.snps.vcf.gz | bcftools norm -m+ | bcftools view --threads 1 -m2 -M2 -Oz -o vcf/${sampleID}.chr${i}.snps.het.vcf.gz
     tabix vcf/${sampleID}.chr${i}.snps.het.vcf.gz) &
    pids[${j}]=$!
    j=$((j+1))
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done

# concatenate individual chromosome VCFs
ls -d -1 "${PWD}"/vcf/${sampleID}*.het.vcf.gz > vcf/${sampleID}.vcf.list
vcf-concat --files vcf/${sampleID}.vcf.list > vcf/${sampleID}.snps.het.vcf

# clean up directory
rm vcf/${sampleID}.chr*.vcf*
rm vcf/${sampleID}.vcf.list

## align RNA-seq data to reference with STAR
# parameters inspired by https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md
STAR \
--runMode alignReads \
--runThreadN ${threads} \
--twopassMode Basic \
--genomeDir ${starIndex} \
--varVCFfile vcf/${sampleID}.snps.het.vcf \
--waspOutputMode SAMtag \
--readFilesIn ${fastqDir}/${batch}/${sampleID}_R1_001.fastq.gz ${fastqDir}/${batch}/${sampleID}_R2_001.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outFileNamePrefix ${workDir}/${sampleID}/${sampleID} \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--outFilterMatchNmin 0 \
--limitSjdbInsertNsj 1200000 \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs None \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMattrRGline ID:${sampleID} SM:${sampleID} \
--outSAMattributes NH HI AS nM NM ch \
--chimOutJunctionFormat 0 \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType Junctions WithinBAM SoftClip \
--chimMainSegmentMultNmax 1 \
--genomeLoad NoSharedMemory
