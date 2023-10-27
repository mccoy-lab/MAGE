library(data.table)
library(tidyverse)
library(pbmcapply)
library(cowplot)
library(khroma)
library(broom)

### sQTLs - intron level ###

cs_per_gene_intron <- fread("~/Dropbox/projects/rna_seq_1kg/sQTL_finemapping.introns.txt.gz")

panel_cs_per_gene_intron <- ggplot(data = cs_per_gene_intron, aes(x = factor(numCredibleSets), fill = "fill")) +
  geom_bar(aes(y = after_stat(count))) +
  theme_bw() +
  xlab("Number of credible causal sets") +
  ylab("Number of sIntrons") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  geom_segment(x = 3, xend = 11, y = 3346, yend = 3346, color = "darkgray", lty = "dashed",
               arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +
  annotate(geom = "text", y = 12000, x = 7, label = "4425 sIntrons (3.9%) with\n≥2 credible sets")

sqtl_cs <- fread("~/Dropbox/projects/rna_seq_1kg/sQTL_finemapping.significantAssociations.txt.gz") %>%
  group_by(., intronID, variantCredibleSet) %>%
  summarize(., numVariants = n()) %>%
  as.data.table()

panel_var_per_cs_intron <- ggplot(data = sqtl_cs, aes(x = numVariants, fill = "fill")) +
  geom_histogram(breaks = seq(0.5, 3902.5, 1)) +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, 4000)) +
  #ylim(0, 4200) +
  geom_hline(yintercept = 7718, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 8000, x = 700, label = "7718 eQTLs with a single\ncredible causal variant")

panel_var_per_cs_log_hist_intron <- ggplot(data = sqtl_cs, aes(x = numVariants, fill = "fill")) +
  geom_histogram(breaks = seq(0.5, 3902.5, 1)) +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  coord_cartesian(xlim = c(1, 4000)) +
  #ylim(0, 4200) +
  geom_hline(yintercept = 7718, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 8000, x = 700, label = "7718 eQTLs with a single\ncredible causal variant")

var_per_cs_dt_intron <- layer_data(panel_var_per_cs_intron) %>% as.data.table()

panel_var_per_cs_log_intron <- ggplot(data = var_per_cs_dt_intron, aes(x = x, y = y)) +
  geom_point(color = "#668A99") +
  geom_line() +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  theme(legend.position = "none") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  coord_cartesian(xlim = c(1, 4000)) +
  #ylim(0, 4400) +
  geom_hline(yintercept = 7718, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 7000, x = 100, label = "7718 intron-level sQTLs with a single\ncredible causal variant")

### sQTLs - gene level ###

cs_per_gene_merged <- fread("~/Dropbox/projects/rna_seq_1kg/sQTL_finemapping.gene.mergedCredibleSets.txt.gz")

panel_cs_per_gene_merged <- ggplot(data = cs_per_gene_merged, aes(x = numMergedCredibleSets, fill = "fill")) +
  geom_bar(aes(y = after_stat(count))) +
  theme_bw() +
  xlab("Number of credible causal sets") +
  ylab("Number of sGenes") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  geom_segment(x = 3, xend = 92, y = 1546, yend = 1546, color = "darkgray", lty = "dashed",
               arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +
  annotate(geom = "text", y = 1900, x = 50, label = "3490 sGenes (45.2%) with\n≥2 credible sets")

sqtl_cs_merged <- fread("~/Dropbox/projects/rna_seq_1kg/sQTL_finemapping.gene.mergedCredibleSets.significantAssociations.txt.gz") %>%
  group_by(., geneID, mergedCredibleSet) %>%
  summarize(., numVariants = n()) %>%
  as.data.table()
 
panel_var_per_cs_merged <- ggplot(data = sqtl_cs_merged, aes(x = numVariants, fill = "fill")) +
  geom_histogram(breaks = seq(0.5, 4329.5, 1)) +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, 4400)) +
  #ylim(0, 4200) +
  geom_hline(yintercept = 3569, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 2000, x = 700, label = "3569 eQTLs with a single\ncredible causal variant")

panel_var_per_cs_log_hist_merged <- ggplot(data = sqtl_cs_merged, aes(x = numVariants, fill = "fill")) +
  geom_histogram(breaks = seq(0.5, 4329.5, 1)) +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  coord_cartesian(xlim = c(1, 4400)) +
  #ylim(0, 4200) +
  geom_hline(yintercept = 3569, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 2000, x = 700, label = "3569 eQTLs with a single\ncredible causal variant")

var_per_cs_dt_merged <- layer_data(panel_var_per_cs_merged) %>% as.data.table()

panel_var_per_cs_log_merged <- ggplot(data = var_per_cs_dt_merged, aes(x = x, y = y)) +
  geom_point(color = "#668A99") +
  geom_line() +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  theme(legend.position = "none") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  coord_cartesian(xlim = c(1, 4400)) +
  #ylim(0, 4400) +
  geom_hline(yintercept = 3569, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 3200, x = 100, label = "3569 gene-level sQTLs with a single\ncredible causal variant")

plot_grid(panel_cs_per_gene_intron, panel_var_per_cs_log_intron, 
          panel_cs_per_gene_merged, panel_var_per_cs_log_merged, 
          nrow = 2, labels = c("A.", "B.", "C.", "D."),
          align = "hv", axis = "lrtb")




