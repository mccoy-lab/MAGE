#!/usr/bin/bash

#==========#
# Get args #
#==========#

kgpVCFdir=$1 # The directory containing 1000 Genomes VCF files (one per chromosome) from the New York Genome Center 8/5/2020 release 
             # with names such as CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz
gtexVCF=$2 # The path to the GTEx VCF file
afgrVCF=$3 # The path to the AFGR MKK VCF file
outDir=$4

mafThreshold=0.05 # MAF threshold to test for associations

#=========================================================================#
# Using 1000 Genomes data, compute LD in windows and remove rare variants #
#=========================================================================#

cd ${outDir}
mkdir 1kgp

for i in {1..22}
do
  plink --vcf ${kgpVCFdir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz \
    --keep keep_sorted.txt \
    --make-bed \
    --snps-only \
    --keep-allele-order \
    --maf ${mafThreshold} \
    --indep-pairwise 200 20 0.2 \
    --out ${outDir}/1kgp/1kgp_chr${i}_maf &
done

wait

for i in {1..22}
do
plink --bfile 1kgp/1kgp_chr${i}_maf \
  --extract 1kgp/1kgp_chr${i}_maf.prune.in \
  --make-bed \
  --out 1kgp/1kgp_chr${i}_maf_extractLDindep
done

rm -f merge_bed_list.txt
# merge chromosomes
for i in {2..22}
do
echo 1kgp/1kgp_chr${i}_maf_extractLDindep
done >> merge_bed_list.txt

plink --bfile 1kgp/1kgp_chr1_maf_extractLDindep \
--make-bed \
--merge-list merge_bed_list.txt \
--out 1kgp_autosomal_maf_extractLDindep

cut -f2 1kgp_autosomal_maf_extractLDindep.bim > 1kgp_snp_list.txt

#===============================================================#
# Convert GTEx VCF to PLINK format and extract overlapping SNPs #
#===============================================================#

plink --vcf ${gtexVCF} \
--make-bed \
--out GTEx

cat GTEx.bim | awk '{print $1"\t"$1":"$4":"$6":"$5"\t"0"\t"$4"\t"$5"\t"$6}' > new.bim
mv GTEx.bim GTEx_prev.bim
mv new.bim GTEx.bim

plink --bfile GTEx \
  --extract 1kgp_snp_list.txt \
  --make-bed \
  --snps-only \
  --out GTEx_1kgp_olap

cut -f2 GTEx_1kgp_olap.bim | sort > gtex_snp_list.txt

#===============================================================#
# Convert AFGR VCF to PLINK format and extract overlapping SNPs #
#===============================================================#

plink --vcf ${afgrVCF} \
--make-bed \
--out AFGR

cat AFGR.bim | awk '{print $1"\t"$1":"$4":"$6":"$5"\t"0"\t"$4"\t"$5"\t"$6}' > new.bim
mv AFGR.bim AFGR_prev.bim
mv new.bim AFGR.bim

plink --bfile AFGR \
  --extract 1kgp_snp_list.txt \
  --make-bed \
  --snps-only \
  --out AFGR_1kgp_olap

cut -f2 AFGR_1kgp_olap.bim | sort > afgr_snp_list.txt

join -1 1 -2 1 gtex_snp_list.txt afgr_snp_list.txt > snp_list.txt

plink --bfile 1kgp_autosomal_maf_extractLDindep \
  --extract snp_list.txt \
  --make-bed \
  --out 1kgp_merge_ready

plink --bfile GTEx_1kgp_olap \
  --extract snp_list.txt \
  --make-bed \
  --snps-only \
  --out gtex_merge_ready

plink --bfile AFGR_1kgp_olap \
  --extract snp_list.txt \
  --make-bed \
  --snps-only \
  --out afgr_merge_ready

#============================#
# Combine the three datasets #
#============================#

echo "gtex_merge_ready" > merge_bed_list.txt
echo "afgr_merge_ready" >> merge_bed_list.txt

plink --bfile 1kgp_merge_ready \
  --merge-list merge_bed_list.txt \
  --make-bed \
  --biallelic-only \
  --out merged_datasets

#============================================#
# Run PCA and ADMIXTURE on the combined data #
#============================================#

plink --bfile merged_datasets \
  --pca

for i in 6 7 8 9 10
do 
admixture ${i} -j8 > admixture_log_${i}.txt &
done
