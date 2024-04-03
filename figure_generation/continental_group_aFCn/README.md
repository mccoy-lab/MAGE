# Continental group `log2(aFC)` comparison

The script in this directory generates the following figures in the MAGE manuscript:

- Fig. S27
- Fig. S28

## Stratified _cis_-eQTL mapping

We applied the MAGE [eQTL mapping pipeline](https://github.com/mccoy-lab/MAGE/tree/main/analysis_pipeline) to each continental group within the dataset, performing _cis_-eQTL mapping, fine-mapping, and effect size estimation with aFCn.

Our QTL mapping scripts were identical to those on the MAGE Github except for the stratification by continental group. Further details are described in the MAGE supplement (section **18.2**).

## Comparison of `log2(aFC)`

The `afcn_comparison.R` script compares log2(aFC) effect sizes calculated within each continental group in MAGE. Further details are described in the MAGE supplement (section **18.2**). There are two main components of the script:

#### 1. Identify shared causal signals between pairs of continental groups.

- We define a shared causal signal as two fine-mapped credible sets (CSs) that uniquely share at least one variant.

- In the process, we remove "multi-matching" credible sets, where the CS in the focal continental group shares variants with >1 CS in the comparison continental group or vice versa.

#### 2. Compare aFCn between shared causal signals.

- We calculate LD between the lead eQTL of each credible set and a variant that is shared between the credible sets of both continental groups. This allows us to reverse the aFCn sign for lead eQTLs that are in negative LD.

- We perform a two-sided Welch's t-test to determine if log2(aFC) is significantly different between the two continental groups.

- We perform Bonferroni correction across all continental groups to account for multiple testing.