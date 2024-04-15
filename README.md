<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/mccoy-lab/MAGE/main/images/MAGE_logo.large_no_bg_white_letters_w_outline.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/mccoy-lab/MAGE/main/images/MAGE_logo.large_no_bg_black_letters_w_outline.png">
  <img alt="MAGE logo" src="https://raw.githubusercontent.com/mccoy-lab/MAGE/main/images/MAGE_logo.large_white_bg_black_letters.png">
</picture>

# MAGE: Multi-ancestry Analysis of Gene Expression

[![DOI](https://zenodo.org/badge/451943672.svg)](https://zenodo.org/doi/10.5281/zenodo.10072080)

MAGE comprises RNA-seq data from lymphoblastoid cell lines derived from 731 individuals from the [1000 Genomes Project (1KGP)](https://doi.org/10.1038/nature15393), representing 26 globally-distributed populations across five continental groups. These data offer a large, geographically diverse, open access resource to facilitate studies of the distribution, genetic underpinnings, and evolution of variation in human transcriptomes and include data from several ancestry groups that were poorly represented in previous studies.


## Data Access

### Raw reads:
Newly generated RNA sequencing data for the 731 individuals (779 total libraries) is available on the Sequence Read Archive (Accession: **[PRJNA851328](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA851328)**).

### Processed data
Processed gene expression matrices and QTL mapping results (as well as a host of other downstream data) are currently available on Zenodo (**[MAGEv1.0 Zenodo link](https://zenodo.org/doi/10.5281/zenodo.10535719)**; **upload in progress**) as well as Dropbox (**[MAGEv1.0 Dropbox link](https://www.dropbox.com/s/yc2ux5v0zoq1sbt/MAGE.v1.0.data.zip?dl=0)**).


Briefly, this repo contains the following data:
1. Sample metadata and sequencing metrics
2. Gene expression and splicing matrices used for e/sQTL mapping and analyses of global trends of expression/splicing diversity
3. cis-e/sQTL mapping results, including aFC estimates for cis-eQTLs
4. Functional annotations of cis-e/sQTLs
5. Results of colocalization analysis between MAGE e/sQTLs and complex trait GWAS from the [PAGE](https://doi.org/10.1038/s41586-019-1310-4) study
6. Results of analyses of global trends of expression/splicing diversity
7. Jointly-generated top genotype PCs for samples in MAGE and other resources with paired WGS/RNA-seq data (Geuvadis, GTEx, AFGR)

READMEs are provided for all data in the repo.

If you're having trouble accessing the data or if there's other data you'd like, feel free to contact us and we'd be happy to share what you need over Globus.

### Variant calls

The high-coverage variant calls used for QTL mapping were previously generated by the [New York Genome Center (NYGC)](https://doi.org/10.1016/j.cell.2022.08.004) and are available through the 1KGP [FTP site](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/).


## Code

Code used for data processing and downstream analyses is made available in the [analysis_pipeline/](analysis_pipeline/) directory, along with READMEs describing how each script is run.

Code used to produce major figures/panels in the manuscript is made available in the [figure_generation/](figure_generation/) directory.


## The MAGE manuscript

For more information about the MAGE resource as well as analyses performed using this resource, please see our paper:

**[Sources of gene expression variation in a globally diverse human cohort](https://www.biorxiv.org/content/10.1101/2023.11.04.565639)**<br>
Dylan J. Taylor, Surya B. Chhetri, Michael G. Tassia, Arjun Biddanda, Stephanie M. Yan, Genevieve L. Wojcik, Alexis Battle, Rajiv C. McCoy†

### Citing MAGE

If you use MAGE data in your own work, please cite the paper linked above.
