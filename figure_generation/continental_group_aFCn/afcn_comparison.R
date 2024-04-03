library(tidyverse)
library(data.table)
library(pbmcapply)
library(khroma)
library(wesanderson)
library(RColorBrewer)

#======================#
# Required input files #
#======================#

# Prefix of `eQTL_finemapping.allAssociations.<contGroup>.txt` file from `run_susie.sh`
finemap_prefix <- "eQTL_finemapping.allAssociations."
# Prefix of `eQTL_finemapping.leadVariant_aFCn.<contGroup>.txt` file from `run_susie.sh`
leadvar_prefix <- "eQTL_finemapping.leadVariant_aFCn.noX."
# Prefix of aFCn results files `aFCn.results.<contGroup>.txt`
afcn_prefix <- "aFCn.results."
# Directory of MAGE continent-specific VCFs
mage_vcfs <- "mage_VCFs/"
# Directory to write VCFs after subsetting to lead/fine-mapped variants
ld_vcfs <- "ld_vcfs/"

# List of continental groups
conts <- c("AFR", "AMR", "EAS", "EUR", "SAS")
# Unique pairs of two continents (order agnostic)
cont_comb <- apply(combn(conts, 2), 2, paste, collapse = "_")
# Unique pairs of two continents (order matters)
cont_perm <- c(cont_comb,
               apply(combn(rev(conts), 2), 2, paste, collapse = "_")) %>%
  sort()
# Number of MAGE samples in each continental group
cont_sizes <- as.list(c(196, 113, 141, 142, 139))
names(cont_sizes) <- conts

# Path to bcftools view install
bcftools <- "/home/syan11/code/bcftools-1.18/bcftools view"
# Path to plink install
plink <- "/home/syan11/code/plink_1.90/plink"


#============================================================#
# Find credible sets that can be compared between continents #
#============================================================#

# Read in all eQTL finemapping sets across continents
all_finemap <- pbmclapply(as.list(conts),
                          function(x) {
                            # read in and remove chrX variants
                            fread(cmd = paste0("grep -v 'X:' ", finemap_prefix, x, ".txt")) %>%
                              # remove results with no credible set
                              .[!(is.na(variantCredibleSet))] %>%
                              # create a unique ID from gene + CS ID (since some genes have multiple)
                              .[, unique_id := paste0(geneID, "_", variantCredibleSet)] %>%
                              # annotate with continent
                              .[, continent := x]
                          }) %>%
  rbindlist()

# Determine whether a given credible set overlaps between two populations
find_cs_overlap <- function(in_id, pop_set, all_finemap) {
  # Separate out continent IDs
  pop1_id <- strsplit(pop_set, "_")[[1]][1]
  pop2_id <- strsplit(pop_set, "_")[[1]][2]
  
  # Get gene from gene-CS ID
  gene <- strsplit(in_id, "_")[[1]][1]
  
  # Variants in credible set from pop1
  pop1 <- all_finemap[unique_id == in_id & continent == pop1_id]
  
  # Get all credible sets for this gene from pop2
  cs_list <- unique(all_finemap[geneID == gene & continent == pop2_id]$variantCredibleSet)
  
  if (length(cs_list) == 0) { # No credible sets in pop2
    out <- data.frame(geneID = gene,
                      pop1_id = in_id,
                      pop2_id = NA,
                      overlap = "No CS in other continent")
    # If there is at least one credible set in pop2, check for overlap
  } else {
    # Merge population dataframes by variant ID
    merged <- merge(pop1, all_finemap[geneID == gene & continent == pop2_id],
                    by = "variantID", all.y = TRUE) %>%
      # Summarize which pop2 CSs sucessfully merged
      group_by(geneID.y, unique_id.x, unique_id.y) %>%
      summarize() %>%
      # Reduce NA lines and keep top row for each group (pop2 ID)
      group_by(geneID.y, unique_id.y) %>%
      summarize(unique_id.x = first(unique_id.x),
                unique_id.y = first(unique_id.y)) %>%
      setDT() %>%
      # Annotate with whether CSs overlap
      .[, overlap := "CSs share variants"] %>%
      .[is.na(unique_id.x), overlap := "CSs do not share variants"] %>%
      .[, unique_id.x := in_id] %>%
      setnames(c("geneID", "pop2_id", "pop1_id", "overlap")) %>%
      .[, c("geneID", "pop1_id", "pop2_id", "overlap")]
    
    # If more than one pop2 credible set overlaps
    if (nrow(merged[overlap == "CSs share variants"]) > 1) {
      out <- data.frame(geneID = gene,
                        pop1_id = in_id,
                        pop2_id = paste0(merged$pop2_id, collapse = ","),
                        overlap = "Multiple CS matches in other continent")
    } else { # Otherwise, write output
      out <- merged
    }
  }
  return(out)
}
# Loop over all pairs of continents and find pairwise overlapping credible sets
cs_overlaps <- pbmclapply(as.list(cont_perm),
                          function(cont_set) {
                            cont1 <- strsplit(cont_set, "_")[[1]][1]
                            cs_ids <- unique(all_finemap[continent == cont1]$unique_id)
                            pbmclapply(as.list(cs_ids),
                                       function(x) find_cs_overlap(x, cont_set, all_finemap)) %>%
                              rbindlist() %>%
                              .[, conts := cont_set]
                          }) %>%
  rbindlist()
fwrite(cs_overlaps, "cs_overlaps.txt", sep = "\t")

# What percent of CSs are comparable between CGs?
cs_overlaps %>%
  group_by(conts) %>%
  summarize(perc_overlap = sum(overlap == "CSs share variants") / n()) %>%
  setorder(perc_overlap)


#=====================================================#
# Identify CSs where other CG has multiple CS matches #
#=====================================================#

# Get vector of CSs that share at least one variant
# Some of these are cases where the CG2 CS matches multiple CG1 CSs,
# but these appear as "CSs share variants" results in CG1 and need to be removed separately
roi <- which(cs_overlaps$overlap == "CSs share variants")

# Get table of CSs where there are multiple matches to another cont
# The matched CSs need to be removed from `roi`
mult_match <- cs_overlaps[overlap == "Multiple CS matches in other continent"] %>%
  # Generate population pair ID for the opposite comparison
  .[, c("cont1", "cont2") := tstrsplit(conts, "_")] %>%
  .[, rev_conts := paste0(cont2, "_", cont1)]
# Get `roi` indices that need to be removed
mult_subset <- pbmclapply(as.list(roi),
                          function(x) {
                            row <- cs_overlaps[x]
                            cont <- row$conts
                            id <- row$pop2_id
                            if (id %in% mult_match$pop1_id) {
                              if (nrow(mult_match[pop1_id == id & rev_conts == cont]) > 0) {
                                return(FALSE)
                              } else {
                                return(TRUE)
                              }
                            } else {
                              return(TRUE)
                            }
                          }) %>%
  unlist()
# Remove multi-matching indices from `roi`
roi <- roi[mult_subset]
fwrite(data.frame("roi" = roi), "roi.txt", col.names = FALSE)

# Annotate multi-matching CSs
cs_overlaps[roi[!mult_subset], overlap := "Multiple CS matches in this continent"]

# Reorganize and relabel data for plotting CS overlaps
cs_overlaps_plot <- cs_overlaps %>%
  group_by(conts, pop1_id) %>%
  arrange(overlap) %>%
  # Collapse into unique CSs per focal continental group
  summarize(olap = paste0(overlap, collapse = ",")) %>%
  setDT() %>%
  .[olap %like% "CSs do not share variants", summary := "Gene has >=1 CS in comparison group; none share variants with focal CS"] %>%
  # If focal CS shares variants with any comparison CS
  .[olap %like% "CSs share variants", summary := "CS shares variants with exactly 1 CS in comparison group"] %>%
  # If comparison CS shares variants with multiple focal CSs (including this one)
  .[olap %like% "Multiple CS matches in this continent", summary := "1 CS in comparison group shares variants with multiple focal CSs"] %>%
  .[olap == "Multiple CS matches in other continent", summary := "CS shares variants with multiple CSs in comparison group"] %>%
  .[olap == "No CS in other continent", summary := "Gene has no CS in comparison group"]
# Reorder and relabel categories
cs_overlaps_plot$summary <- factor(cs_overlaps_plot$summary,
                                   levels = c("Gene has no CS in comparison group",
                                              "Gene has >=1 CS in comparison group; none share variants with focal CS",
                                              "1 CS in comparison group shares variants with multiple focal CSs",
                                              "CS shares variants with multiple CSs in comparison group",
                                              "CS shares variants with exactly 1 CS in comparison group"),
                                   labels = c("Gene has no CS in comparison group",
                                              "Gene has >=1 CS in comparison group;\nnone share variants with focal CS",
                                              "1 CS in comparison group shares variants\nwith multiple focal CSs",
                                              "CS shares variants with\nmultiple CSs in comparison group",
                                              "CS shares variants with\nexactly 1 CS in comparison group"))

# Barplot of CS overlaps between continental group pairs
ggplot(test, aes(x = conts, fill = summary)) +
  geom_bar() +
  scale_fill_manual(values = as.vector(color("light")(5))) +
  scale_x_discrete(labels = gsub("_", "\n", cont_perm)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(margin = margin(10,0,0,0)),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size = 9),
        legend.title = element_blank()) +
  labs(y = "Number of credible sets")
ggsave("cs_overlap.pdf", width = 11, height = 4.8)


#===============================================================#
# Prepare data to compare aFCn between comparable credible sets #
#===============================================================#

# Read in all aFCn results
afcn <- lapply(as.list(conts),
               function(x) {
                 fread(paste0(afcn_prefix, x, ".txt")) %>%
                   .[, cont := x]
               }) %>%
  rbindlist()

# Read in file of all lead variants (to match credible set ID to aFCn lead variant)
lead_vars <- lapply(as.list(conts),
                    function(x) {
                      fread(paste0(leadvar_prefix, x, ".txt")) %>%
                        .[, cont := x]
                    }) %>%
  rbindlist()

# Lookup table of all variant IDs in lead variant VCFs
leadvar_list <- pbmclapply(as.list(conts),
                           function(x) {
                             fread(cmd = paste0("grep -v '##' ", ld_vcfs, x, ".leadVars.vcf")) %>%
                               .[, c("ID")] %>%
                               .[, cont := x]
                           }) %>%
  rbindlist()

# Subset MAGE VCFs to all lead variants, for calculating LD quickly in next section
fwrite(lead_vars[, c("variantID")] %>% unique(), # Variant IDs for bcftools
       paste0("lead_vars.txt"), col.names = FALSE)
# Subset VCFs with bcftools
lapply(as.list(conts),
       function(x) system(command = paste(bcftools,
                                          "-o", paste0(ld_vcfs, x, ".leadVars.vcf"),
                                          "-i 'ID=@lead_vars.txt'",
                                          paste0(mage_vcfs, "MAGE.", x, ".Xdip.vcf.gz"))))

# Subset MAGE VCFs to all finemapped variants, for calculating LD quickly in next section
# (leadVars VCF is smaller and runs faster for LD calculation, but doesn't include all variants)
fwrite(all_finemap[, c("variantID")] %>% unique(), # Variant IDs for bcftools
       paste0("all_finemap_variants.txt"), col.names = FALSE)
# Subset VCFs with bcftools
lapply(as.list(conts),
       function(x) system(command = paste(bcftools,
                                          "-o", paste0(ld_vcfs, x, ".finemapVars.vcf"),
                                          "-i 'ID=@all_finemap_variants.txt'",
                                          paste0(mage_vcfs, "MAGE.", x, ".Xdip.vcf.gz"))))


#=====================#
# Run aFCn comparison #
#=====================#

# Function to calculate LD (r, not r2) between a pair of lead variants
# lead1/2: Lead variant IDs
# pop1: Continental group to calculate LD in
# invcf: VCF to calculate LD with
# ld_prefix: Prefix of plink output
get_plink_r <- function(lead1, lead2, pop1, gene, invcf, ld_prefix) {
  # Set names of files to read in/create
  outprefix <- paste0(ld_prefix, pop1, "_", lead1, "_", lead2, "_", gene) # Prefix for output files
  outvcf <- paste0(outprefix, ".vcf") # Output VCF of two variants of interest
  outld <- paste0(outprefix, ".ld") # plink LD output file
  
  # Subset VCF to just the two variants of interest
  system(command = paste(bcftools, "-o", outvcf,
                         "-i", paste0("'ID==", paste0('"', lead1, '"'),
                                      " || ID==", paste0('"', lead2, '"'), "'"),
                         invcf))
  # Calculate r with plink (using system2 to suppress plink output)
  system2(command = plink,
          args = c("--r square", "--keep-allele-order",
                   paste("--vcf", outvcf),
                   paste("--out", outprefix)),
          stdout = FALSE)
  # Get LD
  ld <- fread(outld)[1, V2]
  # Remove intermediate files
  system(command = paste0("rm ", outprefix, "*"))
  return(ld)
}

# Function to calculate LD between two credible sets (to determine whether they're in negative LD)
# lead1/2: Lead variant IDs
# pop1/2: Continental groups of credible sets
# gene: Gene ID
# all_finemap: `all_finemap` dataframe
# subset: Row of `cs_overlaps`
calculate_ld <- function(lead1, lead2, pop1, pop2, gene, all_finemap, subset) {
  if (lead1 == lead2) { # If lead variants are the same, LD is 1
    p1_ld <- 1
    p2_ld <- 1
  # Otherwise, pick an overlapping variant & calculate LD between lead and overlap SNPs in each continent
  } else {
    p1_vars <- all_finemap[continent == pop1 & unique_id == subset$pop1_id] # All variants in CG1
    p2_vars <- all_finemap[continent == pop2 & unique_id == subset$pop2_id] # All variants in CG2
    # Merge credible sets to find an overlapping variant
    merged <- merge(p1_vars, p2_vars, by = "variantID")
    
    # Handle case where only one variant overlaps
    if (nrow(merged) == 1) {
      if (merged[1, ]$variantID == lead1) { # If lead variant 1 is the one overlap
        p1_ld <- 1 # Set LD in CG1 to 1
        # Get LD between lead1 and lead2 in CG2
        p2_ld <- get_plink_r(lead1, lead2, pop2, gene,
                             paste0(ld_vcfs, pop2, ".leadVars.vcf"), ld_vcfs)
      } else if (merged[1, ]$variantID == lead2) { # If lead variant 2 is the one overlap
        # Get LD between lead1 and lead2 in CG1
        p1_ld <- get_plink_r(lead1, lead2, pop1, gene,
                             paste0(ld_vcfs, pop1, ".leadVars.vcf"), ld_vcfs)
        p2_ld <- 1 # Set LD in CG2 to 1
      } else { # Otherwise, calculate LD normally
        # Get first overlapping variant
        overlap_var <- merged[(variantID != lead1) & (variantID != lead2)][1, ]$variantID
        # LD with `finemapVars` VCF because `leadVars` VCF may not include the overlapping variant
        p1_ld <- get_plink_r(lead1, overlap_var, pop1, gene,
                             paste0(ld_vcfs, pop1, ".finemapVars.vcf"), ld_vcfs)
        p2_ld <- get_plink_r(lead2, overlap_var, pop2, gene,
                             paste0(ld_vcfs, pop2, ".finemapVars.vcf"), ld_vcfs)
      }
    # Handle case where only two variants overlap and one or both may be a lead variant
    } else if (nrow(merged) == 2) {
      # If the only overlapping variants are the two lead variants, calculate LD between them
      if (nrow(merged[(variantID != lead1) & (variantID != lead2)]) == 0) {
        p1_ld <- get_plink_r(lead1, lead2, pop1, gene,
                             paste0(ld_vcfs, pop1, ".leadVars.vcf"), ld_vcfs)
        p2_ld <- get_plink_r(lead1, lead2, pop2, gene,
                             paste0(ld_vcfs, pop2, ".leadVars.vcf"), ld_vcfs)
      # Otherwise, calculate LD but make sure we use the correct variant ID
      } else {
        p1_ld <- get_plink_r(lead1, merged[(variantID != lead1)]$variantID[1], pop1, gene,
                             paste0(ld_vcfs, pop1, ".finemapVars.vcf"), ld_vcfs)
        p2_ld <- get_plink_r(lead2, merged[(variantID != lead2)]$variantID[1], pop2, gene,
                             paste0(ld_vcfs, pop2, ".finemapVars.vcf"), ld_vcfs)
      }
    # If >2 variants overlap, calculate LD between lead and an overlapping variant
    } else {
      # Get ID of the first overlapping variant, excluding lead variants
      overlap_var <- merged[(variantID != lead1) & (variantID != lead2)][1, ]$variantID
      p1_ld <- get_plink_r(lead1, overlap_var, pop1, gene,
                           paste0(ld_vcfs, pop1, ".finemapVars.vcf"), ld_vcfs)
      p2_ld <- get_plink_r(lead2, overlap_var, pop2, gene,
                           paste0(ld_vcfs, pop2, ".finemapVars.vcf"), ld_vcfs)
    }
  }
  return(c(p1_ld, p2_ld))
}

# Function to use Welch's t-test to compare aFCs from two continents
# x1/2: Means of two aFCs; se1/2: standard errors of two aFCs
welch_ttest <- function(x1, x2, se1, se2, pop1, pop2) {
  # Get vi from population size (used to calculate v)
  n1 <- cont_sizes[[pop1]]
  n2 <- cont_sizes[[pop2]]
  v1 <- n1 - 1
  v2 <- n2 - 1
  
  # Calculate t statistic
  t <- (x1 - x2) / sqrt(se1^2 + se2^2)
  
  # Calculate degrees of freedom
  s1 <- se1 * sqrt(n1) # Sample standard deviation
  s2 <- se2 * sqrt(n2)
  v <- (s1^2/n1 + s2^2/n2)^2 / (s1^4/(n1^2*v1) + s2^4/(n2^2*v2))
  
  # Calculate two-tailed p-value (1 - the difference between the one-sided pvalues)
  pval <- 1 - abs(pt(t, v) - pt(-t, v))
}

# Function to compare aFCn intervals for two credible sets
compare_afcn <- function(idx, cs_overlaps, lead_vars, afcn) {
  subset <- cs_overlaps[idx, -c("overlap")] # Row of `cs_overlaps` to pull info from
  # Get gene name
  gene <- subset$geneID
  # Get continent names
  pop1 <- strsplit(subset$conts, "_")[[1]][1]
  pop2 <- strsplit(subset$conts, "_")[[1]][2]
  # Get credible set IDs
  cs1 <- strsplit(subset$pop1_id, "_")[[1]][2]
  cs2 <- strsplit(subset$pop2_id, "_")[[1]][2]
  
  # Get lead variants from `lead_vars` data
  lead1 <- lead_vars[geneID == gene & cont == pop1 & variantCredibleSet == cs1]$variantID
  lead2 <- lead_vars[geneID == gene & cont == pop2 & variantCredibleSet == cs2]$variantID
  # Determine LD between credible sets
  ld_out <- calculate_ld(lead1, lead2, pop1, pop2, gene, all_finemap, subset)
  p1_ld <- ld_out[1]
  p2_ld <- ld_out[2]
  ld_sign <- p1_ld * p2_ld # Negative if lead variants are in negative LD
  
  # Get aFCn results from `afcn` data
  afcn1 <- afcn[gene_id == gene & cont == pop1 & variant_id == lead1]
  afcn2 <- afcn[gene_id == gene & cont == pop2 & variant_id == lead2]
  
  # Get aFC (x) and aFC standard error (se)
  x1 <- afcn1$log2_aFC
  x2 <- afcn2$log2_aFC
  se1 <- afcn1$log2_aFC_error
  se2 <- afcn2$log2_aFC_error
  # Calculate p-value from Welch's t-test
  if (ld_sign < 0) { # If LD is negative, flip the sign of aFC for other continent
    x2 <- -afcn2$log2_aFC
    pval <- welch_ttest(x1, x2, se1, se2, pop1, pop2)
  } else {
    pval <- welch_ttest(x1, x2, se1, se2, pop1, pop2)
  }
  
  # Output original data plus aFCn comparison data
  return(cbind(subset, data.frame(pop1_lead = lead1,
                                  pop1_logafc = x1,
                                  pop2_lead = lead2,
                                  pop2_logafc = x2,
                                  welch_pval = pval,
                                  pop1_ld = p1_ld,
                                  pop2_ld = p2_ld)))
}
# Look for aFCn interval overlaps
afcn_overlaps <- pbmclapply(as.list(roi),
                            function(x) compare_afcn(x, cs_overlaps, lead_vars, afcn)) %>%
  rbindlist()

# Perform Bonferroni correction on Welch's t-test p-values
afcn_overlaps$welch_padj <- p.adjust(afcn_overlaps$welch_pval, method = "bonferroni")
# Annotate with overlap by Welch's t-test
afcn_overlaps$welch_olap <- afcn_overlaps$welch_padj >= 0.05
fwrite(afcn_overlaps, "afcn_overlaps.txt", sep = "\t")


#========================#
# Plotting aFCn overlaps #
#========================#

### Heatmap of aFC overlap percentage

# Summarize overlapping/non-overlapping intervals for heatmap
heatmap_dt <- copy(afcn_overlaps) %>%
  # Only keep unique continent combinations for plotting
  .[conts %in% cont_comb] %>%
  # Make separate columns for each continent
  .[, c("cont1", "cont2") := tstrsplit(conts, "_")] %>%
  .[, c("cont1", "cont2", "welch_olap")] %>%
  # Summarize overlap for each pair of continents
  group_by(cont1, cont2) %>%
  summarize(perc_overlap = round(sum(as.logical(welch_olap), na.rm = TRUE) / n(), digits = 3),
            n_olap = sum(as.logical(welch_olap), na.rm = TRUE),
            n_total = n())
# Reorder y axis so that heatmap is in bottom left corner
heatmap_dt$cont2 <- factor(heatmap_dt$cont2, levels = rev(conts))

ggplot(heatmap_dt, aes(x = cont1, y = cont2, fill = perc_overlap)) +
  geom_tile(color = "black", linetype = 1, lwd = 0.2) +
  geom_text(aes(label = perc_overlap)) +
  theme_classic() +
  scale_fill_gradient2(low = "#d73027", mid = "#fffdbf", high = "#4575b4",
                       midpoint = 0.95,
                       limits = c(0.9, 1),
                       breaks = c(0.9, 1),
                       na.value = NA) +
  labs(fill = "Fraction of\ncredible sets\nwith same\nallelic fold change\n") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        line = element_blank())
ggsave("afcn_overlap_heatmap.pdf", width = 6.4, height = 5)

### Plot of aFCn by aFCn for each pair of continental groups

# Reorder true/false levels
afcn_overlaps$welch_olap <- factor(afcn_overlaps$welch_olap, levels = c(TRUE, FALSE))

# Make cont-cont identity lines so that plot is symmetrical
sample <- afcn_overlaps[1:5,]
sample$conts <- c("AFR_AFR", "AMR_AMR", "EAS_EAS", "EUR_EUR", "SAS_SAS")
afcn_overlaps_plot <- rbind(afcn_overlaps, sample)
afcn_overlaps_plot$conts <- factor(afcn_overlaps_plot$conts,
                                   levels = c("AFR_AFR", "AMR_AFR", "EAS_AFR", "EUR_AFR", "SAS_AFR",
                                              "AFR_AMR", "AMR_AMR", "EAS_AMR", "EUR_AMR", "SAS_AMR",
                                              "AFR_EAS", "AMR_EAS", "EAS_EAS", "EUR_EAS", "SAS_EAS",
                                              "AFR_EUR", "AMR_EUR", "EAS_EUR", "EUR_EUR", "SAS_EUR",
                                              "AFR_SAS", "AMR_SAS", "EAS_SAS", "EUR_SAS", "SAS_SAS"))

ggplot(afcn_overlaps_plot) +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(aes(x = pop1_logafc, y = pop2_logafc),
              method = 'lm', color = "darkgray") +
  geom_point(aes(x = pop1_logafc, y = pop2_logafc, color = welch_olap),
             size = 0.5, alpha = 0.5) +
  scale_color_manual(values = as.vector(color("light")(2))) +
  facet_wrap(~ conts) +
  theme_bw() +
  # theme(strip.background = element_blank(), # remove facet headers for figure
  #       strip.text.x = element_blank()) +
  labs(x = expression(paste(italic(log)[2], "aFC")),
       y = expression(paste(italic(log)[2], "aFC")),
       color = "Same aFC")
ggsave("afc_by_afc.pdf", width = 10, height = 8.8)