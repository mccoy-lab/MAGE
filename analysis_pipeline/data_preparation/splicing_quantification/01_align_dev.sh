#!/bin/bash

module load bcftools
module load htslib
module load vcftools
module load STAR

cd /scratch16/rmccoy22/1KGP_expression/star_alignments

threads=24
our_id=`sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list.txt`
coriell_id=`echo ${our_id} | cut -f1 -d '-'`
participant_id=`cut -f1,3- /data/rmccoy22/1KGP_RNASEQ/updated_metadata.tsv | grep -P "${our_id}\t" | cut -f2`
batch="batch"`cut -f1,3- /data/rmccoy22/1KGP_RNASEQ/updated_metadata.tsv | grep -P "${our_id}\t" | cut -f3`
sex=`cut -f1,3- /data/rmccoy22/1KGP_RNASEQ/updated_metadata.tsv | grep -P "${our_id}\t" | cut -f7`

mkdir /scratch16/rmccoy22/1KGP_expression/star_alignments/output/${our_id}

## generate individual sample-specific VCF files with heterozygous variants
j=1
if [ ${sex} = "female" ]
then
    (i="X"
     echo "Processing chromosome ${i}."
     vcf_file=/data/rmccoy22/1KGP_VCF_3202/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.eagle2-phased.v2.vcf.gz
     bcftools view --threads ${threads} --no-update -s ${participant_id} -v snps ${vcf_file} | bcftools view --threads ${threads} --no-update -e 'GT=".|."' -Oz -o vcf/${our_id}.chr${i}.snps.vcf.gz
     bcftools view --threads ${threads} --no-update -i 'GT="het"' vcf/${our_id}.chr${i}.snps.vcf.gz | bcftools norm -m+ | bcftools view --threads ${threads} -m2 -M2 -Oz -o vcf/${our_id}.chr${i}.snps.het.vcf.gz
     tabix vcf/${our_id}.chr${i}.snps.het.vcf.gz) &
    pids[${j}]=$!
    j=$((j+1))
fi
for i in {1..22}; do
    (echo "Processing chromosome ${i}."
     vcf_file=/data/rmccoy22/1KGP_VCF_3202/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz
     bcftools view --threads 1 --no-update -s ${participant_id} -v snps ${vcf_file} | bcftools view --threads 1 --no-update -e 'GT=".|."' -Oz -o vcf/${our_id}.chr${i}.snps.vcf.gz
     bcftools view --threads 1 --no-update -i 'GT="het"' vcf/${our_id}.chr${i}.snps.vcf.gz | bcftools norm -m+ | bcftools view --threads 1 -m2 -M2 -Oz -o vcf/${our_id}.chr${i}.snps.het.vcf.gz
     tabix vcf/${our_id}.chr${i}.snps.het.vcf.gz) &
    pids[${j}]=$!
    j=$((j+1))
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done

# concatenate individual chromosome VCFs
ls -d -1 "${PWD}"/vcf/${our_id}*.het.vcf.gz > vcf/${our_id}.vcf.list
vcf-concat --files vcf/${our_id}.vcf.list > vcf/${our_id}.snps.het.vcf

# clean up directory
rm vcf/${our_id}.chr*.vcf*
rm vcf/${our_id}.vcf.list

## align RNA-seq data to reference with STAR
# parameters inspired by https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md
STAR \
--runMode alignReads \
--runThreadN ${threads} \
--twopassMode Basic \
--genomeDir /scratch16/rmccoy22/1KGP_expression/star_alignments/ref/star_index_oh149 \
--varVCFfile vcf/${our_id}.snps.het.vcf \
--waspOutputMode SAMtag \
--readFilesIn /data/rmccoy22/1KGP_RNASEQ/raw_reads/${batch}/${our_id}_R1_001.fastq.gz /data/rmccoy22/1KGP_RNASEQ/raw_reads/${batch}/${our_id}_R2_001.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outFileNamePrefix /scratch16/rmccoy22/1KGP_expression/star_alignments/output/${our_id}/${our_id} \
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
--outSAMattrRGline ID:${our_id} SM:${our_id} \
--outSAMattributes NH HI AS nM NM ch \
--chimOutJunctionFormat 0 \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType Junctions WithinBAM SoftClip \
--chimMainSegmentMultNmax 1 \
--genomeLoad NoSharedMemory

chgrp -R rmccoy22 /scratch16/rmccoy22/1KGP_expression/star_alignments/output/${our_id}
