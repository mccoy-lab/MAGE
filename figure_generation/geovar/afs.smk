#!python3

import numpy as np
import pandas as pd

import pickle, gzip
from tqdm import tqdm
from pathlib import Path
from io import StringIO

# ---- Parameters for inference in MAGE data ---- #
vcf_file = "/scratch16/rmccoy22/1KGP_expression/sharedData/sample_vcfs/1KGP_731-samples_all.filtered.phased.vcf.gz"
nominal_eqtls = "/scratch16/rmccoy22/1KGP_expression/QTL_analysis/eQTL_discovery/eQTL_mapping_results/fastQTL_significant_results/fastqtl_all.significantPairs.txt.gz"
finemapped_eqtls = "/scratch16/rmccoy22/1KGP_expression/QTL_analysis/eQTL_discovery/eQTL_mapping_results/susie_results/eQTL_finemapping.significantAssociations.txt.gz"
finemapped_sqtls = "/scratch16/rmccoy22/1KGP_expression/QTL_analysis/sQTL_discovery/sQTL_mapping_results/susie_results/sQTL_finemapping.gene.mergedCredibleSets.bestHits.txt.gz"


# ------- Rules Section ------- #
localrules:
    all,


rule all:
    input:
        "results/af_tables/nominal_mage.eqtl.tsv",
        "results/af_tables/finemapped_mage.eqtl.tsv",
        "results/af_tables/finemapped_single_mage.sqtl.tsv",
        "results/af_tables/finemapped_single_mage.eqtl.tsv",
        "results/variant_lists/gtex_caviar_finemapped_eqtl_single.txt.gz",


rule filter_finemapping_cred_sets:
    """Filter credible sets to only their single most significant variant per-gene."""
    input:
        finemapped_eqtls=finemapped_eqtls,
    output:
        variants="results/variant_lists/finemapped_eqtl_single.txt.gz",
    run:
        df = pd.read_csv(input.finemapped_eqtls, sep="\s")
        tot_df = (
            df.groupby(["geneID", "variantCredibleSet"])[["variantPIP", "variantID"]]
            .agg(np.max)
            .reset_index()[["geneID", "variantID"]]
        )
        tot_df.to_csv(output.variants, sep="\t", index=None)


rule assign_allele_freq_nominal:
    """Assign allele frequencies to the nominal eQTLs"""
    input:
        vcf_file=vcf_file,
        nominal_eqtls=nominal_eqtls,
    output:
        snplist=temp("variants.nominal.mage_eqtl.txt"),
        af_tsv="results/af_tables/nominal_mage.eqtl.tsv",
    shell:
        """
        zcat {input.nominal_eqtls} | awk '{{print $2}}' | sort | uniq > {output.snplist}
        bcftools query -i \"ID=@{output.snplist}\" -f \"%CHROM\t%ID\t%REF\t%ALT\t%AC\t%AF\t%AF_AFR\t%AF_EUR\t%AF_SAS\t%AF_EAS\t%AF_AMR\n\" {input.vcf_file} | awk \'BEGIN{{print \"CHR\tSNP\tA1\tA2\tMAC\tMAF\tAFR\tEUR\tSAS\tEAS\tAMR\"}}; {{print}};\' > {output.af_tsv}
        """


rule assign_allele_freq_finemapped_total:
    """Assign allele frequencies to the nominal eQTLs"""
    input:
        vcf_file=vcf_file,
        finemapped_eqtls=finemapped_eqtls,
    output:
        snplist=temp("variants.finemapped.mage_eqtl.txt"),
        af_tsv="results/af_tables/finemapped_mage.eqtl.tsv",
    shell:
        """
        zcat {input.finemapped_eqtls} | awk '{{print $2}}' | sort | uniq > {output.snplist}
        bcftools query -i \"ID=@{output.snplist}\" -f \"%CHROM\t%ID\t%REF\t%ALT\t%AC\t%AF\t%AF_AFR\t%AF_EUR\t%AF_SAS\t%AF_EAS\t%AF_AMR\n\" {input.vcf_file} | awk \'BEGIN{{print \"CHR\tSNP\tA1\tA2\tMAC\tMAF\tAFR\tEUR\tSAS\tEAS\tAMR\"}}; {{print}};\' > {output.af_tsv}
        """


rule assign_allele_freq_finemapped_single_total:
    """Assign allele frequencies to the nominal eQTLs"""
    input:
        vcf_file=vcf_file,
        finemapped_single_eqtls="results/variant_lists/finemapped_eqtl_single.txt.gz",
    output:
        snplist=temp("variants.finemapped.single.mage_eqtl.txt"),
        af_tsv="results/af_tables/finemapped_single_mage.eqtl.tsv",
    shell:
        """
        zcat {input.finemapped_single_eqtls} | awk '{{print $2}}' | sort | uniq > {output.snplist}
        bcftools query -i \"ID=@{output.snplist}\" -f \"%CHROM\t%ID\t%REF\t%ALT\t%AC\t%AF\t%AF_AFR\t%AF_EUR\t%AF_SAS\t%AF_EAS\t%AF_AMR\n\" {input.vcf_file} | awk \'BEGIN{{print \"CHR\tSNP\tA1\tA2\tMAC\tMAF\tAFR\tEUR\tSAS\tEAS\tAMR\"}}; {{print}};\' > {output.af_tsv}
        """


rule assign_allele_freq_top_sqtl:
    input:
        vcf_file=vcf_file,
        sqtls=finemapped_sqtls,
    output:
        snplist=temp("variants.finemapped.single.mage_sqtl.txt"),
        af_tsv="results/af_tables/finemapped_single_mage.sqtl.tsv",
    shell:
        """
        zcat {input.sqtls} | awk '{{print $2}}' | sort | uniq > {output.snplist}
        bcftools query -i \"ID=@{output.snplist}\" -f \"%CHROM\t%ID\t%REF\t%ALT\t%AC\t%AF\t%AF_AFR\t%AF_EUR\t%AF_SAS\t%AF_EAS\t%AF_AMR\n\" {input.vcf_file} | awk \'BEGIN{{print \"CHR\tSNP\tA1\tA2\tMAC\tMAF\tAFR\tEUR\tSAS\tEAS\tAMR\"}}; {{print}};\' > {output.af_tsv}
        """


rule isolate_gtex_eqls:
    input:
        "/data/rmccoy22/GTEx_Analysis_v8_eQTL/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF_with_Allele.txt.gz",
    output:
        "results/variant_lists/gtex_caviar_finemapped_eqtl_single.txt.gz",
    shell:
        """
        zcat {input} | awk \'NR==1 {{print \"CHROM\tPOS\tREF\tALT\tPIP\tTISSUE\tGENE\"}}; NR > 1 {{split($3, x, \"_\"); OFS=\"\t\"; print x[1], x[2], x[3], x[4], $6, $1,$2}}\' | gzip > {output} 
        """


rule isolate_kgpex_dapg_replicates:
    """Isolate KGPEX variants that replicate and do not replicate in DAPG sets """
    input:
        vcf_file=vcf_file,
        non_replicates="results_forDylan/mage_dapg_metaCS_nonreplicates_nonzeroPIP.top_eVariants.090123.tsv",
        replicates="results_forDylan/mage_dapg_metaCS_replicates_nonzeroPIP.top_eVariants.090123.tsv",
    output:
        snplist_replicates="results/variant_lists/mage_dapg_metaCS_replicates_snplist.txt",
        snplist_nonreplicates="results/variant_lists/mage_dapg_metaCS_nonreplicates_snplist.txt",
        af_replicates_tsv="results/af_tables/mage_dapg_metaCS_replicates.top_eVariants.tsv",
        af_nonreplicates_tsv="results/af_tables/mage_dapg_metaCS_nonreplicates.top_eVariants.tsv",
    shell:
        """
        awk \' NR > 1 {{print $2}}\' {input.replicates} | sort | uniq > {output.snplist_replicates}
        awk \' NR > 1 {{print $2}}\' {input.non_replicates} | sort | uniq > {output.snplist_nonreplicates}
        bcftools query -i \"ID=@{output.snplist_replicates}\" -f \"%CHROM\t%ID\t%REF\t%ALT\t%AC\t%AF\t%AF_AFR\t%AF_EUR\t%AF_SAS\t%AF_EAS\t%AF_AMR\n\" {input.vcf_file} | awk \'BEGIN{{print \"CHR\tSNP\tA1\tA2\tMAC\tMAF\tAFR\tEUR\tSAS\tEAS\tAMR\"}}; {{print}};\' > {output.af_replicates_tsv}
        bcftools query -i \"ID=@{output.snplist_nonreplicates}\" -f \"%CHROM\t%ID\t%REF\t%ALT\t%AC\t%AF\t%AF_AFR\t%AF_EUR\t%AF_SAS\t%AF_EAS\t%AF_AMR\n\" {input.vcf_file} | awk \'BEGIN{{print \"CHR\tSNP\tA1\tA2\tMAC\tMAF\tAFR\tEUR\tSAS\tEAS\tAMR\"}}; {{print}};\' > {output.af_nonreplicates_tsv}
        """
