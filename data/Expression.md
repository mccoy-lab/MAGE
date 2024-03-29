# Expression and global trend data for MAGE v0.1

# Table of contents

* [**Table of contents**](#table-of-contents)
* [**Overview**](#overview)
* [**Expression**](#expression)
	* [Gene expression (counts)](#gene-expression-counts)
	* [Gene expression (mu)](#gene-expression-mu)
	* [Gene expression (vst)](#gene-expression-vst)
* [**Expression variation**](#expression-variation)
    * [Proportion of variance explained by population](#proportion-of-variance-explained-by-population)
    * [Population permutation test](#population-permutation-test)
    * [Proportion of variance explained by sequencing batch](#proportion-of-variance-explained-by-batch)
* [**Relationship between fixation index and differential expression**](#Relationship-between-fixation-index-and-differential-expression)

<br><br>

# Overview

<img src="/images/dropbox.png" width="15" style="float: bottom;"> **[MAGEv.01 Expression and global trend data]()**

Gene expression matrices and global trends in gene expression results for the MAGE v0.1 data set are publicly available for download through Dropbox.

<br><br>

# Expression

<img src="/images/dropbox.png" width="15" style="float: bottom;"> **[MAGEv.01 Expression matrices]()**

Gene expression quantification was performed using [Salmon](https://salmon.readthedocs.io/en/latest/) for genes in [GENCODE v38](https://www.gencodegenes.org/human/release_38.html). Transcript-level estimates was converted into gene-level estimates using [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html). Technical replicates for the 24 samples sequenced in triplicate were collapsed using the [`collapseReplicates`](https://rdrr.io/bioc/DESeq2/man/collapseReplicates.html) function included in DESeq2.

<br><br>

## Gene expression (counts)

`expression.counts.csv`: Matrix containing the number of reads mapped (obtained from `NumReads` column in each samples [`quant.sf`](https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats) output from `Salmon`) for each gene (rows) and each sample (columns).

<br><br>

## Gene expression (mu)

`expression.mu.csv`: Modelled mean expression for each gene (rows) and sample (column) generated by [`DESeq2`](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseq2-model). Mean expression incorporates both sample- and gene-level normalization factors.

<br><br>

## Gene expression (VST)

`expression.vst.csv`: Blind variance stabilized transformation (VST) gene-level count data generated using the [`vst`](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variance-stabilizing-transformation) command in `DESeq2`. As the VST matrix is normalized across the expression matrix, accounting for normalization factors, gene-wise dispersion, and removed dependence of variance on the mean without *a priori* knowledge of the design matrix (`blind=TRUE`), these data were used in subsequent analyses to estimate global trends in gene expression. Columns correspond to samples and rows to gene.

<br><br>

# Expression variation

<img src="/images/dropbox.png" width="15" style="float: bottom;"> **[MAGEv.01 Expression variation]()**

Using the VST expression matrix described above, we quantified the expression trends at the biological scale of populations as well as the technical scale of sequence batches. 

## Proportion of variance explained by population

Two-stage regression results for estimating the proportion of variance explained by *continentalGroup* label (e.g., SAS & AMR) and *population* label (e.g., LWK & JPT). After regressing out the effects of *sex* and *batch* on expression (using the VST expression matrix), residuals were then regressed, separately, against *continentalGroup* and *population*. Proportion of variation explain was then calculated as the SSR/SST.

In summary, the following linear models are used for regression:
<br>
$$VST(expression) \sim sex + batch + u$$
$$u \sim continentalGroup + v_1$$
$$u \sim population + v_2$$

Where *u* are the residuals from the first-stage regression, and *v* are the residuals from the second-stage regressions.

A summary of these results are presented in `PVE_regressions.csv`. The columns in this file are as follows:
1. `gene`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID
2. `stage1.SSR_sex`: Regression sum of squares for *sex* fixed-effect in stage-one regression
3. `stage1.SSR_batch`: Regression sum of squares for *batch* fixed-effect in stage-one regression
4. `stage1.SSR_res`: Residual sum of squares ($u$) in stage-one regression
5. `stage2a.SSR_continentalGroup`: Regression sum of squares for *continentalGroup* fixed-effect in stage-two regression
6. `stage2a.SSR_res`: Residual sum of squares ($v_1$) in stage-two regression 
7. `stage2b.SSR_population`: Regression sum of squares for *population* in stage-two regression
8. `stage2b.SSR_res`: Residual sum of squares ($v_2$) in stage-two regression
9. `continentalGroupPVE`: Proportion of variation explained by *continentalGroup* fixed-effect ($`SSR_{continentalGroup}/SST`$)
10. `populationPVE`: Proportion of variation explained by *continentalGroup* fixed-effect ($`SSR_{population}/SST`$)

<br><br>

## Population permutation test

To test whether the proportion of variance in gene expression explained by continental group and/or population was greater than expected by chance, we leveraged permutation tests by shuffling *continentalGroup* and *population* labels (with replacement) and repeated the [two-stage ANOVA above](#proportion-of-variance-explained-by-population) to compute a null distribution of expected proportion of variance explained (*n = 1000 replicates*).

`PVE_permutations.csv` contains a summary of these permutation results, where the mean of proportion of variance explained across all genes was computed for each permutation. Each row of this table corresponds to an independent permutation procedure, and the two columns are:
1. `mean_continentalGroupPVE`: Mean proportion of variance explained by *continentalGroup* label across all genes (e.g., the mean of column (*9.*) in `PVE_regressions.csv`).
2. `mean_populationPVE`: Mean proportion of variance explained across all genes by *population* label (e.g., the mean of column (*10.*) in `PVE_regressions.csv`).

<br><br>

## Proportion of variance explained by sequencing batch

We estimated the proportion of variation explained by sequencing *batch* and *sample* using the 24 samples sequenced in triplicate. For this test, we used a VST matrix of autosomal gene counts generated from the raw count matrix without running [`collapseReplicates`](#expression). A type-2 ANOVA was performed using the following linear-regression model:
<br>
$$VST(expression) \sim batch + sample$$  

The proportion of variance explained by *batch* and *sample* were computed as the regression sum of squares (for *batch* or *sample*, respectively) divided by the total sum of squares.

`Replicate_ANOVA.csv` summarises the results from these tests. The columns of this file are:
1. `Gene`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID
2. `SSRbatch`: Regression sum of squares for *batch* categorical fixed-effect
3. `SSRsample`: Regression sum of squares for *sample* categorical fixed-effect
4. `SSresiduals`: Residual sum of squares
5. `ANOVAp_batch`: Type-2 ANOVA p-value for *batch* variable
6. `ANOVAp_sample`: Type-2 ANOVA p-value for *sample* variable
7. `qvalue_batch`: q-value for *batch* variable (FDR = 0.05)
8. `qvalue_sample`: q-value for *sample* variable (FDR = 0.05)

<br><br>

# Relationship between fixation index and differential expression

<img src="/images/dropbox.png" width="15" style="float: bottom;"> **[MAGEv.01 Fst vs DE]()**

To investigate the relationship between differential gene expression and population differentiation in allele frequency, we computed Weir & Cockerham's $F_{st}$ using [`vcflib`](https://github.com/vcflib/vcflib/blob/master/doc/wcFst.md) for each finemapped eQTL. These results are presented in `Finemapped_Fst.csv` with the following columns:
1. `geneID`: The [Ensembl](https://useast.ensembl.org/Homo_sapiens/Info/Index) gene ID
2. `continentalGroup`: Foreground population for $F_{st}$ calculation
3. `mean_wcFst`: Mean $F_{st}$ across all finemapped eQTLs for the focal eGene, where the foreground allele frequency is measured for the population measured in the `continentalGroup` column, and the background allele frequency is measured in the remaining four `continentalGroup` categories (e.g., AFR vs. EUR + SAS + EAS + AMR). 
4. `DEpadj`: Differential gene expression adjusted p-value (FDR = 0.05) computed for focal population against mean expression of all other populations (e.g., AFR vs. EUR + SAS + EAS + AMR)
