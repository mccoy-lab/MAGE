# code to process data for PCA
# to compare ancestry composition
# of various datasets

cd /scratch16/rmccoy22/rmccoy22/kgpex_pca
mkdir 1kgp

# compute LD in windows and remove rare variants

for i in {1..22}
do
  plink --vcf ~/data_rmccoy22/1KGP_VCF_3202/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz \
    --keep keep_sorted.txt \
    --make-bed \
    --snps-only \
    --keep-allele-order \
    --maf 0.05 \
    --indep-pairwise 200 20 0.2 \
    --out 1kgp/1kgp_chr${i}_maf5percent &
done

wait

# extract LD-independent SNPs
for i in {1..22}
do
plink --bfile 1kgp/1kgp_chr${i}_maf5percent \
  --extract 1kgp/1kgp_chr${i}_maf5percent.prune.in \
  --make-bed \
  --out 1kgp/1kgp_chr${i}_maf5percent_extractLDindep
done

rm -f merge_bed_list.txt
# merge chromosomes
for i in {2..22}
do
echo 1kgp/1kgp_chr${i}_maf5percent_extractLDindep
done >> merge_bed_list.txt

plink --bfile 1kgp/1kgp_chr1_maf5percent_extractLDindep \
--make-bed \
--merge-list merge_bed_list.txt \
--out 1kgp_autosomal_maf5percent_extractLDindep

cut -f2 1kgp_autosomal_maf5percent_extractLDindep.bim > 1kgp_snp_list.txt

# convert GTEx to PLINK file and extract overlapping SNPs

nohup plink --vcf phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz \
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

# convert MKK AGFR data to PLINK file and extract overlapping SNPs

plink --vcf mkk/MKK.filtered.biallelicAndRare.autosomal.chr.id.vcf.gz \
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

### extract overlapping SNPs

join -1 1 -2 1 gtex_snp_list.txt afgr_snp_list.txt > snp_list.txt

plink --bfile 1kgp_autosomal_maf5percent_extractLDindep \
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

# combine data

echo "gtex_merge_ready" > merge_bed_list.txt
echo "afgr_merge_ready" >> merge_bed_list.txt

plink --bfile 1kgp_merge_ready \
  --merge-list merge_bed_list.txt \
  --make-bed \
  --biallelic-only \
  --out merged_datasets

# run PCA and ADMIXTURE on the combined data
plink --bfile merged_datasets \
  --pca

for i in 6 7 8 9 10
do 
nohup ./dist/admixture_linux-1.3.0/admixture merged_datasets.bed ${i} -j8 > admixture_log_${i}.txt &
done
