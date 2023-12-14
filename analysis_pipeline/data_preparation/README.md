# MAGE Data Preparation

This directory contains all the components of the data preparation pipeline used as part of the initial MAGE publication. This pipeline is broadly split into three sections:
1. Preparing MAGE VCFs from 1KGP VCFs
2. Quantifying expression from raw reads using Salmon
3. Alignment with STAR and splicing quantification with Leafcutter

## Preparing MAGE VCFs

Code to produce MAGE VCFs is in the [`subset_VCFs/`](subset_VCFs/) directory.
