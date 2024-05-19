# Quantifying expression variance

The `R` script contained in this directory (`generate_expression_ANOVAs.R`) will compute the following expression ANOVAs from the VST expression data:
1. Variance explained by continental group and population labels (`PVE_regressions.csv`)
2. Residual variance in gene expression after regressing out the effects of batch and sex - partitioned by continental group label (`contGroup_variance.csv`)
3. Residual variance in gene expression after regressing out the effects of batch and sex - partitioned by population label (`pop_variance.csv`)

The script uses the following data as input:

- `expression.vst.csv`: [VST transformed expression](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variance-stabilizing-transformation) generated after [collapsing technical replicates](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#collapsing-technical-replicates) with the `DESeq2` package.
- `sample.metadata.MAGE.v1.0.txt`: Sample metadata
- `expression_filteredGenes.MAGE.v1.0.txt.gz`: List of ensembl IDs corresponding to genes with ≥6 counts and ≥0.1 TPM in at least 20% (147/731) of samples.

All input files and copies of the output from `generate_expression_ANOVAs.R` can be found in the [Zenodo Repository](https://zenodo.org/records/10535719).

---

*Authors note:*

*The `R` script used to generate these data is inefficient, as it is not parallelized. As such, we recommend reconstructing the ANOVAs to leverage parallelization such as those available through the `pbmclapply` package.*