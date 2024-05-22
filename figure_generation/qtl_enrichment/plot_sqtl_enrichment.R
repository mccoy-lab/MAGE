#!/usr/bin/env R

# Generate enrichment plots for the sQTLs
#############################################
#                                           #
# ENRICHMENT PLOTS                          #
#                                           #
#############################################


#' Generate Enrichment Plots for sQTLs
#'
#' This function reads statistical data from specified files within a given directory,
#' reshapes it for plotting, and then generates two plots: one showing log2 fold enrichment
#' and another showing the proportion of sQTLs. These plots are saved to a PDF file.
#'
#' @param results_dir A string that specifies the path to the directory containing
#'        the input data files. This directory should include TSV files with binomial
#'        and Fisher's test statistics for sQTLs.
#'
#' @return None explicitly, but this function saves a PDF file containing the plots
#'         in the working directory.
#'
#' @examples
#' generate_enrichment_plots("/path/to/results_dir")
#'
#' @importFrom dplyr filter select left_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_tile scale_fill_manual geom_point coord_flip theme_bw labs geom_hline theme
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob gpar
#' @importFrom data.table fread
#' @export
generate_enrichment_plots <- function(results_dir) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)

  # Read data
  df_binom <- fread(file.path(results_dir, "QTL_1X.vep.loftee.vcf.allAnnot.results.BinomialSTATS.tsv"))
  df_fisher <- fread(file.path(results_dir, "QTL_1X.vep.loftee.vcf.allAnnot.results.fishersSTATS.updated.tsv"))

  # Reshape data for easy plotting
  df_long <- df_fisher %>%
    select(Annotation, log2_Fold_Enrichment = "log2(FoldEnrich)", log10_pvalue = "-log10(P-Value)", Proportion_Observed) %>%
    pivot_longer(cols = -Annotation, names_to = "Variable", values_to = "Value")

  # Split data for different plots
  df_enrichment <- filter(df_long, Variable == "log2_Fold_Enrichment")
  df_proportion <- filter(df_long, Variable == "Proportion_Observed")
  df_significance <- filter(df_long, Variable == "log10_pvalue")

  # Adjust factor levels
  factor_order <- factor(df_enrichment$Annotation, levels = unique(df_enrichment$Annotation))
  df_enrichment$Annotation <- factor(df_enrichment$Annotation, levels = factor_order)
  df_proportion$Annotation <- factor(df_proportion$Annotation, levels = factor_order)

  # Generate enrichment plot
  plot_enrichment <- ggplot(df_enrichment, aes(x = Annotation, y = Value)) +
    geom_tile(aes(fill = rep(c(TRUE, FALSE), length.out = nrow(df_enrichment))), width = 0.9, height = Inf, alpha = 0.25) +
    scale_fill_manual(values = c("grey", "transparent")) +
    geom_point(shape = 23, size = 3, fill = "orange") +
    coord_flip() +
    theme_bw() +
    labs(y = "log2(Fold Enrichment)", x = "") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  # Join significance data and generate proportion plot
  df_proportion <- df_proportion %>%
    left_join(select(df_significance, Annotation, Significance = Value), by = "Annotation") %>%
    mutate(Asterisk = case_when(
      Significance <= 3 ~ "",
      Significance > 3 & Significance < 10 ~ "*",
      Significance >= 10 ~ "**"
    ))

  plot_proportion <- ggplot(df_proportion, aes(x = Annotation, y = Value)) +
    geom_bar(stat = "identity", fill = "orange") +
    geom_text(aes(label = Asterisk, y = Value + 0.15), vjust = 1, hjust = 1.0, size = 5) +
    coord_flip() +
    theme_bw() +
    labs(y = "Proportion of sQTLs", x = "") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_blank())

  # Output to PDF
  pdf("./FoldEnrich_proportion.pdf", width = 8, height = 7)
  grid.arrange(
    arrangeGrob(plot_enrichment, top = textGrob("A", x = 0.05, y = 0.85, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 18))),
    arrangeGrob(plot_proportion, top = textGrob("B", x = 0.05, y = 0.85, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 18))),
    ncol = 2, widths = c(3, 1)
  )
  dev.off()
}


########################################
#                                      #
# WITH CONFIDENCE INTERVALS            #
#                                      #
########################################


#' Generate Enrichment Plots for sQTLs with Confidence Intervals
#'
#' This function reads Fisher's test statistics from a specified file within a given directory,
#' reshapes it for plotting, and then generates two plots: one showing log2 fold enrichment with
#' confidence intervals, and another showing the proportion of sQTLs. These plots are saved to a PDF file.
#'
#' @param results_dir A string specifying the path to the directory containing
#'        the input data files. This directory should include a TSV file with Fisher's test
#'        statistics for sQTLs.
#'
#' @return None explicitly, but this function saves a PDF file containing the plots
#'         in the specified directory.
#'
#' @examples
#' generate_enrichment_plots_with_ci("/path/to/results_dir")
#'
#' @importFrom dplyr filter select left_join mutate
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_replace_all
#' @importFrom ggplot2 ggplot geom_tile scale_fill_manual geom_point geom_errorbar coord_flip theme_bw labs geom_hline theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob gpar
#' @importFrom data.table fread
#' @export
generate_enrichment_plots_with_ci <- function(results_dir) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
  library(grid)

  # Read Fisher's test statistics
  df_fisher <- fread(file.path(results_dir, "QTL_1X.vep.loftee.vcf.allAnnot.results.fishersSTATS.updated.tsv"))

  # Reshape the data for easy plotting
  df_long <- df_fisher %>%
    select(Annotation, log2_Fold_Enrichment = "log2(FoldEnrich)", 
           CI_upper = "CI Upper (Fold Enrichment, log2)",
           CI_lower = "CI Lower (Fold Enrichment, log2)",
           log10_pvalue = "-log10(P-Value)", Proportion_Observed) %>%
    pivot_longer(cols = -Annotation, names_to = "Variable", values_to = "Value")

  # Modify the Annotation column
  df_long$Annotation <- df_long$Annotation %>%
    str_replace_all("Utr", "UTR") %>%
    str_replace_all("Tf", "TF") %>%
    str_replace_all("5Th", "5th") %>%
    str_replace_all("Nmd", "NMD") %>%
    str_replace_all("Lof", "LoF")

  # Split and reshape data for plots
  df_enrichment <- filter(df_long, Variable == "log2_Fold_Enrichment")
  df_proportion <- filter(df_long, Variable == "Proportion_Observed")
  df_significance <- filter(df_long, Variable == "log10_pvalue")
  factor_order <- factor(df_enrichment$Annotation, levels = unique(df_enrichment$Annotation))

  # Extract CI information and join with enrichment data
  df_ci <- df_long %>%
    filter(Variable %in% c("CI_upper", "CI_lower")) %>%
    pivot_wider(names_from = Variable, values_from = Value)
  df_enrichment <- df_enrichment %>%
    left_join(df_ci, by = "Annotation") %>%
    mutate(alt_id = rep(c(TRUE, FALSE), ceiling(nrow(.) / 2))[1:nrow(.)]) %>%
    arrange(factor_order) %>%
    mutate(Annotation = factor(Annotation, levels = factor_order))

  # Create the enrichment plot with error bars
  plot_enrichment <- ggplot(df_enrichment, aes(x = Annotation, y = Value)) +
    geom_tile(aes(fill = alt_id), width = 0.9, height = Inf, alpha = 0.25) +
    scale_fill_manual(values = c("grey", "transparent")) +
    geom_point(shape = 23, size = 3, fill = "orange") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    coord_flip() +
    theme_bw() +
    labs(y = "log2(Fold Enrichment)", x = "") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

  # Generate proportion plot with significance markers
  df_proportion <- df_proportion %>%
    left_join(select(df_significance, Annotation, Significance = Value), by = "Annotation") %>%
    mutate(Asterisk = case_when(
      Significance <= 3 ~ "",
      Significance > 3 & Significance < 10 ~ "*",
      Significance >= 10 ~ "**"),
      Annotation = factor(Annotation, levels = factor_order))

  plot_proportion <- ggplot(df_proportion, aes(x = Annotation, y = Value)) +
    geom_bar(stat = "identity", fill = "orange") +
    geom_text(aes(label = Asterisk, y = Value + 0.15), vjust = 1, hjust = 1.0, size = 5) +
    coord_flip() +
    theme_bw() +
    labs(y = "Proportion of sQTLs", x = "") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_blank())

  # Output to PDF
  pdf(file.path(results_dir, "Plot_FoldEnrich_proportion_with_CI.pdf"), width = 8, height = 7)
  grid.arrange(
    arrangeGrob(plot_enrichment, top = textGrob("", x = 0.05, y = 0.85, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 18))),
    arrangeGrob(plot_proportion, top = textGrob("", x = 0.05, y = 0.85, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 18))),
    ncol = 2, widths = c(3, 1)
  )
  dev.off()
}


########################################
#                                      #
# WITH LOG ODDS RATIO                  #
#                                      #
########################################


#' Generate Log(Odds Ratio) Plots for sQTLs
#'
#' This function reads Fisher's test statistics from a specified file within a given directory,
#' reshapes it for plotting, and generates two types of plots: one showing the log2 of the odds ratio,
#' and another showing the proportion of sQTLs. The results are saved in PDF format.
#'
#' @param results_dir A string specifying the path to the directory containing
#'        the input data files. This directory should include TSV files with Fisher's test
#'        statistics for sQTLs.
#'
#' @return None explicitly, but this function saves a PDF file containing the plots
#'         in the working directory.
#'
#' @examples
#' generate_log_odds_ratio_plots("/path/to/results_dir")
#'
#' @importFrom dplyr filter select left_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_tile scale_fill_manual geom_point coord_flip theme_bw labs geom_hline theme
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob gpar
#' @importFrom data.table fread
#' @export
generate_log_odds_ratio_plots <- function(results_dir) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)

  # Read Fisher's test statistics
  df_fisher <- fread(file.path(results_dir, "QTL_1X.vep.loftee.vcf.allAnnot.results.fishersSTATS.tsv"))

  # Reshape the data for easy plotting
  df_long <- df_fisher %>%
    select(Annotation, log2_Odds_ratio = "log2(Odds-ratio)", log10_pvalue = "-log10(P-Value)", Proportion_Observed) %>%
    pivot_longer(cols = -Annotation, names_to = "Variable", values_to = "Value")

  # Split data for different plots
  df_enrichment <- filter(df_long, Variable == "log2_Odds_ratio")
  df_proportion <- filter(df_long, Variable == "Proportion_Observed")
  df_significance <- filter(df_long, Variable == "log10_pvalue")
  factor_order <- factor(df_enrichment$Annotation, levels = unique(df_enrichment$Annotation))

  # Adjust data tables
  df_enrichment <- df_enrichment[order(factor_order), ]
  df_enrichment$alt_id <- rep(c(TRUE, FALSE), ceiling(nrow(df_enrichment) / 2))[1:nrow(df_enrichment)]
  df_enrichment$Annotation <- factor(df_enrichment$Annotation, levels = factor_order)
  df_proportion$Annotation <- factor(df_proportion$Annotation, levels = factor_order)

  # Create the enrichment plot
  plot_enrichment <- ggplot(df_enrichment, aes(x = Annotation, y = Value)) +
    geom_tile(aes(fill = alt_id), width = 0.9, height = Inf, alpha = 0.25) +
    scale_fill_manual(values = c("grey", "transparent")) +
    geom_point(shape = 23, size = 3, fill = "#FF3333") +
    coord_flip() +
    theme_bw() +
    labs(y = "log2(Odds-ratio)", x = "") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

  # Generate proportion plot with significance markers
  df_proportion <- df_proportion %>%
    left_join(select(df_significance, Annotation, Significance = Value), by = "Annotation") %>%
    mutate(Asterisk = case_when(
      Significance <= 3 ~ "",
      Significance > 3 & Significance < 10 ~ "*",
      Significance >= 10 ~ "**"))

  plot_proportion <- ggplot(df_proportion, aes(x = Annotation, y = Value)) +
    geom_bar(stat = "identity", fill = "#FF3333") +
    geom_text(aes(label = Asterisk, y = Value + 0.15), vjust = 1, hjust = 1.0, size = 5) +
    coord_flip() +
    theme_bw() +
    labs(y = "Proportion of sQTLs", x = "") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_blank())

  # Save output to PDF
  pdf("./OddsRatio_proportion.pdf", width = 8, height = 7)
  grid.arrange(
    arrangeGrob(plot_enrichment, top = textGrob("A", x = 0.05, y = 0.85, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 18))),
    arrangeGrob(plot_proportion, top = textGrob("B", x = 0.05, y = 0.85, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 18))),
    ncol = 2, widths = c(3, 1)
  )
  dev.off()
}


########################################
#                                      #
# GENERATE SCATTER PLOTS               #
#                                      #
########################################


#' Generate Scatter Plot for Fold Enrichment versus P-Value with Bonferroni Correction
#'
#' This function reads Fisher's test statistics from a specified file within a given directory,
#' reshapes it for plotting, and generates an annotated scatter plot displaying log2 fold enrichment
#' against -log10(p-value) with a Bonferroni correction line. The result is saved in PDF format.
#'
#' @param results_dir A string specifying the path to the directory containing
#'        the input data files. This directory should include a TSV file with Fisher's test
#'        statistics for sQTLs.
#' @param num_tests The number of tests performed for Bonferroni correction calculation.
#' @param significance_level The significance level used for the Bonferroni correction, default is 0.05.
#'
#' @return None explicitly, but this function saves a PDF file containing the scatter plot
#'         in the specified directory.
#'
#' @examples
#' generate_scatter_plot("/path/to/results_dir", 26, 0.05)
#'
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_point geom_label_repel geom_hline geom_vline scale_x_continuous theme_bw labs theme
#' @importFrom ggrepel geom_label_repel
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob gpar
#' @importFrom data.table fread
#' @export
generate_scatter_plot <- function(results_dir, num_tests, significance_level = 0.05) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
  library(grid)

  # Read Fisher's test statistics
  df_fisher <- fread(file.path(results_dir, "QTL_1X.vep.loftee.vcf.allAnnot.results.fishersSTATS.updated.tsv"))

  # Reshape the data for easy plotting
  df_fisher <- df_fisher %>%
    select(Annotation, log2_Fold_Enrichment = "log2(FoldEnrich)", 
           CI_upper = "CI Upper (Fold Enrichment, log2)",
           CI_lower = "CI Lower (Fold Enrichment, log2)",
           log10_pvalue = "-log10(P-Value)", Proportion_Observed)

  # Modify the Annotation column
  df_fisher$Annotation <- df_fisher$Annotation %>%
    str_replace_all("Utr", "UTR") %>%
    str_replace_all("Tf", "TF") %>%
    str_replace_all("5Th", "5th") %>%
    str_replace_all("Nmd", "NMD") %>%
    str_replace_all("Lof", "LoF")

  # Calculate Bonferroni corrected p-value threshold
  bonferroni_threshold <- -log10(significance_level / num_tests)

  # Create the scatter plot
  elegant_scatter <- ggplot(df_fisher, aes(x = log2_Fold_Enrichment, y = log10_pvalue)) +
    geom_point(color = "#FF5733", alpha = 0.9) +
    geom_label_repel(aes(label = Annotation), box.padding = 0.35, point.padding = 0.5,
                     segment.color = 'grey50', size = 2.5, fontface = 'bold', fill = 'white') +
    geom_hline(yintercept = bonferroni_threshold, linetype = "dashed", color = "orange") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "orange") +
    scale_x_continuous(breaks = seq(floor(min(df_fisher$log2_Fold_Enrichment)), 
                                    ceiling(max(df_fisher$log2_Fold_Enrichment)), by = 1)) +
    theme_bw() +
    labs(x = "Log2(Fold Enrichment)", y = "-Log10(P-Value)")

  # Output to PDF
  pdf(file.path(results_dir, "Plot_Scatter_FoldEnrich_pvals_bf.pdf"), width = 9, height = 9)
  grid.arrange(
    arrangeGrob(elegant_scatter, top = textGrob("A", x = 0.05, y = 0.85, just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 18))),
    ncol = 1
  )
  dev.off()
}


# Call enrichment plot functions
generate_enrichment_plots("./results_dir")
generate_enrichment_plots_with_ci("./results_dir")
generate_log_odds_ratio_plots("./results_dir")
generate_scatter_plot("./results_dir", 26, 0.05)

