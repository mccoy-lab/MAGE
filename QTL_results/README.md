# eQTL and sQTL results for MAGE v0.1

## Overview

This directory contains eQTL and sQTL mapping results for the MAGE v0.1 data set. For eQTL and sQTL mapping, a nominal pass was performed with [FastQTL](https://github.com/francois-a/fastqtl). To identify eGenes and sGenes, as well as significant nominal associations, a permutation pass was performed, also with [FastQTL](https://github.com/francois-a/fastqtl).

After identifying eGenes and sGenes with FastQTL, fine-mapping was performed with [SuSiE](https://github.com/stephenslab/susieR) to identify independent credible causal sets of variants for each e/sGene.

For MAGE v0.1, QTL mapping is limited to the autosomes; there are no results for chrX or chrY. QTL mapping results for the sex chromosomes will be added in a future release.

All QTL summary statistics are available in tabixed format in the **[MAGE dropbox]()**. All genomic coordinates are relative to the GRCh38 reference. 

## [eQTL results]()

### Filtered genes

A list of the genes that passed expression filtering and were used for eQTL mapping is available in [`expression_filteredGenes.MAGE.v0.1.txt.gz`](), along with it's corresponding tabix index.

The genomic coordinates listed in this file are BED-style: 0-based, half-open.

### eQTL mapping summary

A summary of the eQTL mapping results for the filtered genes is available in [`eQTL_summary.MAGE.v0.1.txt.gz`](). The columns in this file are as follows:
1. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID
2. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol
3. `variantsTested`: The total number of variants tested for an association with expression of the gene
4. `significantAssocs`: The number of variant-gene associations that reached significance in the nominal FastQTL pass
5. `is_eGene`: Whether the gene reached significance in the FastQTL permutation pass (5% FDR)
6. `finemappedCredibleSets`: For the genes that reached significance in the permutation pass (i.e., eGenes), the number of independent credible causal sets identified by fine-mapping with SuSiE

### FastQTL nominal results

The **[eQTL_nominal_results]()** directory contains nominal eQTL results: summary statistics from the nominal and permutation passes of FastQTL. There are three files in this directory:

#### Permutation pass results
Permutation pass results are available in [`eQTL_FastQTL_results.permutation_pass.MAGE.v0.1.txt.gz`](). This file is organized on a by-gene basis: there is one line for each gene in the filtered data set. The columns in this file are as follows:
1. `geneChrom`: The chromosome of the gene
2. `geneChromTSS`: The 1-based position (on `geneChrom`) of the transcription start site (TSS) of the gene
3. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID
4. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol
5. `num_var`: The total number of variants tested for an association with expression of the gene
6. `beta_shape1`: The first parameter value of the fitted beta distribution
7. `beta_shape2`: The second parameter value of the fitted beta distribution
8. `true_df`: The effective degrees of freedom the beta distribution approximation
9. `pval_true_df`: The empirical _p_-value for the beta distribution approximation
10. `top_kgpID`: The ID of the top *nominal* variant, as defined in the [NYGC 1KGP VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
11. `top_rsID`: The dbNSP rsID of the top variant
12. `tss_distance`: The distance from the top variant to the gene transcription start site (TSS)
13. `ma_samples`: The number of samples carrying the minor allele of the top variant
14. `ma_count`: The total number of minor alleles of the top variant across individuals 
15. `maf`: The minor allele frequency of the top variant in the MAGE v0.1 cohort
16. `ref_factor`: A flag indicating if the alternative allele of the top variant is the minor allele in the cohort (1 if AF <= 0.5, -1 if not)
17. `pval_nominal`: The nominal _p_-value of the top variant from linear regression
18. `slope`: The slope of the linear regression withe top variant
19. `slope_se`: The standard error of the slope
20. `pval_perm`: The permutation _p_-value directly obtained from the permutations with the direct method
21. `pval_beta`: The permutation _p_-value obtained via beta approximation
22.  `qval`: Storey q-value derived from `pval_beta` (FDR adjusted). Genes with `qval <= 0.05` are defined as eGenes
23.  `pval_nominal_threshold`: The nominal _p_-value threshold for calling a variant-gene pair significant for the gene

A tabix index is also provided for this file.

#### All nominal pass results
Nominal pass results for ALL tested associations (significant and non-significant) are available in [`eQTL_FastQTL_results.nominal_pass.allAssociations.MAGE.v0.1.txt.gz`](). This file is organized on a by-association basis: there is one line for each tested variant-gene association in the filtered data set. The columns in this file are as follows:
1. `variantChrom`: The chromosome of the tested variant
2. `variantPosition`: The 1-based position (on `variantChrom`) of the tested variant
3. `variantRef`: The reference allele of the tested variant
4. `variantAlt`: The alternative allele of the tested variant
5. `variant_kgpID`: The ID of the tested variant, as defined in the [NYGC 1KGP VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
6. `variant_rsID`: The dbNSP rsID of the tested variant
7. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID of the tested gene
8. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol of the tested gene
9. `tss_distance`: The distance from the tested variant to the TSS of the tested gene
10. `ma_samples`: The number of samples carrying the minor allele of the tested variant
11. `ma_count`: The total number of minor alleles of the tested variant across individuals 
12. `maf`: The minor allele frequency of the tested variant in the MAGE v0.1 cohort
13. `pval_nominal`: The nominal _p_-value of the linear regression
14. `slope`: The slope of the linear regression
15. `slope_se`: The standard error of the slope

A tabix index is also provided for this file.

#### Significant nominal pass results
Significant nominal associations were determined by comparing the nominal _p_-value of the association (recorded in the `pval_nominal` field in the nominal pass results file) to the nominal _p_-value threshold of the corresponding gene determined during the permutation pass (recorded in the `pval_nominal_threshold` field in the permutation pass results file).

Associations for which `pval_nominal < pval_nominal_threshold` were defined as significant nominal associations.

Nominal pass results for ONLY the significant nominal associations are available in [`eQTL_FastQTL_results.nominal_pass.significantAssociations.MAGE.v0.1.txt.gz`](). The columns of this file are as described in the [All nominal pass results]() section above.

A tabix index is also provided for this file.

### SuSiE fine-mapping results

The **[eQTL_finemapping_results]()** directory contains fine-mapping eQTL results: credible sets and summary statistics from SuSiE, and effect sizes from aFCn. Fine-mapping was performed for all genes identified as eGenes in the FastQTL permutation pass. There are four files files in this directory:

#### Credible set summary
Summary statistics for all of the eQTL credible sets identified by SuSiE are available in [`eQTL_finemapping.credibleSet_summary.MAGE.v0.1.txt.gz`](). This file is organized on a by-credible set basis: there is one line for each significant (coverage > 0.95) credible causal set. The columns in this file are as follows:
1. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID of the gene
2. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol of the gene
3. `credibleSet`: The ID of the credible set (within that gene)
4. `numVariants`: The number of variants in that credible set
5. `coverage`: The coverage of the credible set. Essentially, this describes the probability that the credible set contains a true causal variant. We only discovered credible sets with `coverage > 0.95`
6. `minCorr`: The minimum correlation (in genotype) between variants in the credible set
7. `meanCorr`: The mean correlation (in genotype) between variants in the credible set
8. `medianCorr`: The median correlation (in genotype) between variants in the credible set

#### All fine-mapping results
Fine-mapping results for ALL tested variants (included in credible sets and not) are available in [`eQTL_finemapping.allAssociations.MAGE.v0.1.txt.gz`](). This file is organized on a by-association basis: there is one line for each tested variant-gene association in the filtered data set. The columns in this file are as follows:
1. `variantChrom`: The chromosome of the tested variant
2. `variantPosition`: The 1-based position (on `variantChrom`) of the tested variant
3. `variantRef`: The reference allele of the tested variant
4. `variantAlt`: The alternative allele of the tested variant
5. `variant_kgpID`: The ID of the tested variant, as defined in the [NYGC 1KGP VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
6. `variant_rsID`: The dbNSP rsID of the tested variant
7. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID of the tested gene
8. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol of the tested gene
9. `variantPIP`: The posterior inclusion probability (PIP) that the tested variant is causal for the tested gene
10. `variantCredibleSet`: If appropriate, the ID of the credible set to which the tested variant belongs (within the tested gene)

A tabix index is also provided for this file.

#### Significant fine-mapping results
Fine-mapping results for ONLY those variants that were part of a credible causal set are available in [`eQTL_finemapping.significantAssociations.MAGE.v0.1.txt.gz`](). This file is organized as described in the [All fine-mapping results]() section above.

A tabix index is also provided for this file.

#### Lead variant effect sizes
For each credible set, we selected a single "lead" variant with the highest PIP within that credible set. For each lead variant, we calculated its effect on the associated gene as measured by allelic fold change (aFC). The lead variants for each credible set, along with their corresponding aFC effect size is available in [`eQTL_finemapping.leadVariant_aFCn.MAGE.v0.1.txt.gz`](). This file is organized on a by-credible set basis: there is one line for each credible causal set. The columns in this file are as follows:
1. `variantChrom`: The chromosome of the lead variant
2. `variantPosition`: The 1-based position (on `variantChrom`) of the lead variant
3. `variantRef`: The reference allele of the lead variant
4. `variantAlt`: The alternative allele of the lead variant
5. `variant_kgpID`: The ID of the lead variant, as defined in the [NYGC 1KGP VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
6. `variant_rsID`: The dbNSP rsID of the tested variant
7. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID of the tested gene
8. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol of the tested gene
9. `variantPIP`: The posterior inclusion probability (PIP) that the lead variant is causal for the tested gene
10. `variantCredibleSet`: The ID of the credible set to which the lead variant belongs (within the tested gene)
11. `log2_aFC`: An estimate of the effect of the lead variant on expression of the tested gene as measured by log2(aFC). This value is used as the effect size of the lead variant and it's credible set in downstream analyses
12. `log2_aFC_error`: The standard error of the log2(aFC) estimate
13. `log2_aFC_c0`: The estimate of the "reference" expression of the tested gene. This describes the expected expression of the gene (as log2[counts]) on a haplotype that contains the reference allele for causal variants
14. `log2_aFC_min_95_interv`: The lower limit of the 95% CI around the `log2_aFC` estimate
15. `log2_aFC_max_95_interv`: The upper limit of the 95% CI around the `log2_aFC` estimate
16. `log2_aFC_c0_min_95_interv`: The lower limit of the 95% CI around the `log2_aFC_c0` estimate
17. `log2_aFC_c0_max_95_interv`: The upper limit of the 95% CI around the `log2_aFC_c0` estimate

A tabix index is also provided for this file.

## [sQTL results]()

sQTL mapping was done using intron excision ratios from [Leafcutter](https://davidaknowles.github.io/leafcutter/), and was done in two broad steps.

In the first step, we used FastQTL to identify sGenes and significant nominal associations. In a permutation pass, intron excision ratios were regressed onto sample genotypes, using FastQTL's `--phenotype_groups` option to compute a gene-level empirical p-value over all introns of the gene. Then, a nominal pass was run to discover significant variant-intron associations.

In the second step, we used SuSiE to perform fine-mapping on all sGenes identified in the FastQTL permutation pass. We ran SuSiE for each intron of each sGene, to identify independent credible causal sets for each intron. Then, for each sGene, we created gene-level merged credible sets by iteratively merging any overlapping (i.e. contains the same variant) intron-level credible sets for that gene.

### Filtered genes

A list of the genes that passed splicing filtering and were used for sQTL mapping is available in [`splicing_filteredGenes.MAGE.v0.1.txt.gz`](), along with it's corresponding tabix index.

The genomic coordinates listed in this file are BED-style: 0-based, half-open.

### sQTL mapping summary

A summary of the sQTL mapping results for the filtered genes is available in [`sQTL_summary.MAGE.v0.1.txt.gz`](). The columns in this file are as follows:
1. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID
2. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol
3. `splicingClusters`: The number of Leafcutter splicing clusters (after filtering) in the gene
4. `splicingIntrons`: The number of introns (after filtering) across splicing clusters of the gene
5. `variantsTested`: The total number of variants tested for an association with alternative splicing of the gene
6. `significantVariantGeneAssocs`: The number of unique variant-gene associations that reached significance in the nominal FastQTL pass
7. `significantVariantIntronAssocs`: The number of variant-intron associations that reached significance in the nominal FastQTL pass
8. `is_sGene`: Whether the gene reached significance in the FastQTL permutation pass (5% FDR)
9. `intronFinemappedCredibleSets`: For the genes that reached significance in the permutation pass (i.e., sGenes), the sum of the number of fine-mapped credible causal sets for each intron of that gene
10. `mergedFinemappedCredibleSets`: For the genes that reached significance in the permutation pass (i.e., sGenes), the number of merged credible causal sets

### FastQTL nominal results

The **[sQTL_nominal_results]()** directory contains nominal sQTL results: summary statistics from the nominal and permutation passes of FastQTL. There are three files in this directory:

#### Permutation pass results
Permutation pass results are available in [`sQTL_FastQTL_results.permutation_pass.MAGE.v0.1.txt.gz`](). This file largely matches the format of the eQTL permutation pass file described [above](), but for completeness, we describe it here as well.

This file is organized on a by-gene basis: there is one line for each gene in the filtered data set. The columns in this file are as follows:
1. `geneChrom`: The chromosome of the gene
2. `geneChromTSS`: The 1-based position (on `geneChrom`) of the transcription start site (TSS) of the gene
3. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID
4. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol
5. `top_intronID`: The ID of the intron in the top *nominal* association. Intron IDs have the following fields:
	1. The chromosome of the intron
	2. The 1-based end position of the exon upstream of the intron (on the plus strand)
	3. The 1-based start position of the exon downstream of the intron (on the plus strand)
	4. The Leafcutter cluster ID
	5. The Ensembl ID of the gene
6. `num_var`: The total number of variants tested for an association with splicing of the gene
7. `beta_shape1`: The first parameter value of the fitted beta distribution
8. `beta_shape2`: The second parameter value of the fitted beta distribution
9. `true_df`: The effective degrees of freedom the beta distribution approximation
10. `pval_true_df`: The empirical _p_-value for the beta distribution approximation
11. `top_kgpID`: The ID the variant in the top *nominal* association, as defined in the [NYGC 1KGP VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
12. `top_rsID`: The dbNSP rsID of the top variant
13. `tss_distance`: The distance from the top variant to the gene transcription start site (TSS)
14. `ma_samples`: The number of samples carrying the minor allele of the top variant
15. `ma_count`: The total number of minor alleles of the top variant across individuals 
16. `maf`: The minor allele frequency of the top variant in the MAGE v0.1 cohort
17. `ref_factor`: A flag indicating if the alternative allele of the top variant is the minor allele in the cohort (1 if AF <= 0.5, -1 if not)
18. `pval_nominal`: The nominal _p_-value of the top variant from linear regression
19. `slope`: The slope of the linear regression withe top variant
20. `slope_se`: The standard error of the slope
21. `pval_perm`: The permutation _p_-value directly obtained from the permutations with the direct method
22. `pval_beta`: The permutation _p_-value obtained via beta approximation
23.  `group_id`: The ID of the factor used to group tests with the FastQTL `--group_phenotypes` option. This is identical to `ensemblID`
24.  `group_size`: The number of separate phenotypes (introns) that comprise the group (gene)
25.  `qval`: Storey q-value derived from `pval_beta` (FDR adjusted). Genes with `qval <= 0.05` are defined as sGenes
26.  `pval_nominal_threshold`: The nominal _p_-value threshold for calling a variant-intron pair significant for the gene

A tabix index is also provided for this file.

#### All nominal pass results
Nominal pass results for ALL tested associations (significant and non-significant) are available in [`sQTL_FastQTL_results.nominal_pass.allAssociations.MAGE.v0.1.txt.gz`](). As before, this file largely matches the format of the eQTL nominal  pass file described [above](), but for completeness, we describe it here as well.

This file is organized on a by-association basis: there is one line for each tested variant-intron association in the filtered data set. The columns in this file are as follows:
1. `variantChrom`: The chromosome of the tested variant
2. `variantPosition`: The 1-based position (on `variantChrom`) of the tested variant
3. `variantRef`: The reference allele of the tested variant
4. `variantAlt`: The alternative allele of the tested variant
5. `variant_kgpID`: The ID of the tested variant, as defined in the [NYGC 1KGP VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
6. `variant_rsID`: The dbNSP rsID of the tested variant
7. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID of the tested gene
8. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol of the tested gene
9.  `intronID`: The ID of the tested intron. The format of these IDs is described in the [Permutation pass results]() section above
10. `tss_distance`: The distance from the tested variant to the TSS of the tested gene
11. `ma_samples`: The number of samples carrying the minor allele of the tested variant
12. `ma_count`: The total number of minor alleles of the tested variant across individuals 
13. `maf`: The minor allele frequency of the tested variant in the MAGE v0.1 cohort
14. `pval_nominal`: The nominal _p_-value of the linear regression
15. `slope`: The slope of the linear regression
16. `slope_se`: The standard error of the slope

A tabix index is also provided for this file.

#### Significant nominal pass results
Significant nominal associations were determined by comparing the nominal _p_-value of the association (recorded in the `pval_nominal` field in the nominal pass results file) to the nominal _p_-value threshold of the corresponding gene determined during the permutation pass (recorded in the `pval_nominal_threshold` field in the permutation pass results file).

Associations for which `pval_nominal < pval_nominal_threshold` were defined as significant nominal associations.

Nominal pass results for ONLY the significant nominal associations are available in [`sQTL_FastQTL_results.nominal_pass.significantAssociations.MAGE.v0.1.txt.gz`](). The columns of this file are as described in the [All nominal pass results]() section above.

A tabix index is also provided for this file.

### SuSiE fine-mapping results

The **[sQTL_finemapping_results]()** directory contains fine-mapping sQTL results: credible sets and summary statistics from SuSiE. Fine-mapping was performed for all introns of all genes identified as sGenes in the FastQTL permutation pass. There are four files files in this directory:

#### Credible set summary
Summary statistics for all of the sQTL credible sets identified by SuSiE are available in [`sQTL_finemapping.introns.credibleSet_summary.MAGE.v0.1.txt.gz`](). This file is organized on a by-credible set basis: there is one line for each significant (coverage > 0.95) credible causal set. The columns in this file are as follows:
1. `intronID`: The ID of the tested intron. The format of these IDs is described in the [Permutation pass results]() section above
2. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID of the gene
3. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol of the gene
4. `credibleSet`: The ID of the credible set (within that intron)
5. `numVariants`: The number of variants in that credible set
6. `coverage`: The coverage of the credible set. Essentially, this describes the probability that the credible set contains a true causal variant. We only discovered credible sets with `coverage > 0.95`
7. `minCorr`: The minimum correlation (in genotype) between variants in the credible set
8. `meanCorr`: The mean correlation (in genotype) between variants in the credible set
9. `medianCorr`: The median correlation (in genotype) between variants in the credible set

#### All fine-mapping results
Fine-mapping results for ALL tested variants (included in credible sets and not) are available in [`sQTL_finemapping.introns.allAssociations.MAGE.v0.1.txt.gz`](). This file is organized on a by-association basis: there is one line for each tested variant-gene association in the filtered data set. The columns in this file are as follows:
1. `variantChrom`: The chromosome of the tested variant
2. `variantPosition`: The 1-based position (on `variantChrom`) of the tested variant
3. `variantRef`: The reference allele of the tested variant
4. `variantAlt`: The alternative allele of the tested variant
5. `variant_kgpID`: The ID of the tested variant, as defined in the [NYGC 1KGP VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
6. `variant_rsID`: The dbNSP rsID of the tested variant
7. `intronID`: The ID of the tested intron. The format of these IDs is described in the [Permutation pass results]() section above
8. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID of the tested gene
9. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol of the tested gene
10. `variantPIP`: The posterior inclusion probability (PIP) that the tested variant is causal for the tested intron
11. `variantCredibleSet`: If appropriate, the ID of the credible set to which the tested variant belongs (within the tested intron)

A tabix index is also provided for this file.

#### Significant intron-level fine-mapping results
Fine-mapping results for ONLY those variants that were part of a credible causal set are available in [`sQTL_finemapping.introns.significantAssociations.MAGE.v0.1.txt.gz`](). This file is organized as described in the [All fine-mapping results]() section above.

A tabix index is also provided for this file.

#### Significant gene-level fine-mapping results
for each sGene, we created gene-level merged credible sets by iteratively merging any overlapping (i.e. contains the same variant) intron-level credible sets for that gene. These results are available in [`sQTL_finemapping.geneMerged.significantAssociations.MAGE.v0.1.txt.gz`](). This file is organized on a by-association basis: there is one line for each tested variant-gene association in the filtered data set. The columns in this file are as follows:
1. `variantChrom`: The chromosome of the tested variant
2. `variantPosition`: The 1-based position (on `variantChrom`) of the tested variant
3. `variantRef`: The reference allele of the tested variant
4. `variantAlt`: The alternative allele of the tested variant
5. `variant_kgpID`: The ID of the tested variant, as defined in the [NYGC 1KGP VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
6. `variant_rsID`: The dbNSP rsID of the tested variant
7. `ensemblID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID of the tested gene
8. `geneSymbol`: The [HGNC](https://www.genenames.org/) gene symbol of the tested gene
9. `mergedCredibleSet`: The ID of the gene-level merged credible set to which the tested variant belongs (within the tested gene)
10. `intronIDs`: A comma-separated list of all introns that were merged to form the gene-level credible set identified in the `mergedCredibleSet` column

A tabix index is also provided for this file.<br><br>
