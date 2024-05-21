# Analysis of splicing variation within and between human groups

This directory contains code to replicate the splicing analysis ANOVA analyses done for the initial MAGE publication. We used [MANTA](https://github.com/dgarrimar/manta) to perform MANOVA using filtered splicing ratios from Leafcutter. The procedure comprises three broad steps:
1. Run MANTA to quantify the proportion of splicing variance explained by continental group or population for each splicing cluster
2. Run MANTA with continental group and population label permutations (so as to determine significance of proportion of variance distributions)
3. Quantify splicing variance within each continental group and population

*Note: these steps need not be run in any particular order*<br><br>

## Initial splicing ANOVA

The `01_run_manta.sh` script will run MANTA for each of the splicing clusters in the filtered Leafcutter output. For each cluster, we quantify the proportion of splicing variance explained by continental group or population, after controlling for the effects of batch and sex.

The script takes three arguments:
1. `ratiosFile`: Filtered intron excision ratios. This is the `<prefix>_perind.counts.filtered.gz` file created in the [`splicing_quantification`](../data_preparation/splicing_quantification/) step
2. `metadataFile`: A five-column tab-separated file with sample metadata. File should have a header, columns must be: `sampleID`, `batch`, `sex`, `continentalGroup`, and `populations`
3. `outFile`: File to write ANOVA results to

*Note that this script uses 16 threads by default, you are free to change this however you like*

### Output

The output file will have one row for each splicing cluster. The proportion of variance explained by population (after correcting for batch and sex) is `population_SS/batch_sex_residual_SS`. Similarly, the proportion of variance explained by continental group is `continentalGroup_SS/batch_sex_residual_SS`<br><br>

## MANTA permutations

The `02_run_manta_permutations.sh` script will repeat the analysis done by `01_run_manta.sh`, after permuting population and continental group labels, to determine the signficance of the proportions determined from the `01_run_manta.sh` script.

The script takes four arguments:
1. `ratiosFile`: Filtered intron excision ratios. This is the `<prefix>_perind.counts.filtered.gz` file created in the [`splicing_quantification`](../data_preparation/splicing_quantification/) step
2. `metadataFile`: A five-column tab-separated file with sample metadata. File should have a header, columns must be: `sampleID`, `batch`, `sex`, `continentalGroup`, and `populations`
3. `outDir`: Directory to write permutation ANOVA results to
4. `outPrefix`: Prefix of output files. One output file will be created for each permutation.

*Note that the script runs 1000 permutations across 16 threads by default, but you are welcome to change either of these parameters*

### Output

The script will produce one output file for each permutation, matching the output file from the `01_run_manta.sh` script<br><br>

## Splicing variance within groups

The `03_get_per_pop_variance.sh` script will calculate splicing variance for each splicing cluster within 1) each continental group and 2) each population.

The script takes four arguments:
1. `ratiosFile`: Filtered intron excision ratios. This is the `<prefix>_perind.counts.filtered.gz` file created in the [`splicing_quantification`](../data_preparation/splicing_quantification/) step
2. `metadataFile`: A five-column tab-separated file with sample metadata. File should have a header, columns must be: `sampleID`, `batch`, `sex`, `continentalGroup`, and `populations`
3. `outContGroupFile`: File with within-continental group variance quantifications
4. `outPopFile`: File with within-population variance quantifications

### Output

Output files contain splicing variance for each cluster within each group. See the methods section in our manuscript for an explanation of how variance is calculated<br><br>
