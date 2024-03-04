library(data.table)
library(tidyverse)
library(pbmcapply)
library(cowplot)
library(khroma)
library(broom)

### distribution of credible sets per gene ###

cs_per_gene <- fread("~/Downloads/eQTL_finemapping.genes.txt.gz")

panel_cs_per_gene <- ggplot(data = cs_per_gene, aes(x = factor(numCredibleSets), fill = "fill")) +
  geom_bar(aes(y = after_stat(count))) +
  theme_bw() +
  xlab("Number of credible causal sets") +
  ylab("Number of eGenes") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  geom_segment(x = 3, xend = 11, y = 3200, yend = 3200, color = "darkgray", lty = "dashed",
               arrow = arrow(ends = "both", type = "closed", length = unit(0.05, "inches"))) +
  annotate(geom = "text", y = 3200, x = 7, label = "3951 eGenes (26.3%) with\n≥2 credible sets")

### distribution of variants per credible causal set ###

eqtl_cs <- fread("~/Downloads/eQTL_finemapping.credibleSets.txt.gz")

panel_var_per_cs <- ggplot(data = eqtl_cs, aes(x = numVariants, fill = "fill")) +
  geom_histogram(breaks = seq(0.5, 3763.5, 1)) +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, 4000)) +
  ylim(0, 4200) +
  geom_hline(yintercept = 3992, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 4000, x = 700, label = "3992 eQTLs with a single\ncredible causal variant")

panel_var_per_cs_log_hist <- ggplot(data = eqtl_cs, aes(x = numVariants, fill = "fill")) +
  geom_histogram(breaks = seq(0.5, 3763.5, 1)) +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  scale_fill_manual(values = "#668A99") +
  theme(legend.position = "none") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  coord_cartesian(xlim = c(1, 4000)) +
  ylim(0, 4200) +
  geom_hline(yintercept = 3992, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 4000, x = 700, label = "3992 eQTLs with a single\ncredible causal variant")

var_per_cs_dt <- layer_data(panel_var_per_cs) %>% as.data.table()

panel_var_per_cs_log <- ggplot(data = var_per_cs_dt, aes(x = x, y = y)) +
  geom_point(color = "#668A99") +
  geom_line() +
  theme_bw() +
  xlab("Number of variants") +
  ylab("Number of credible sets") +
  theme(legend.position = "none") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  coord_cartesian(xlim = c(1, 4000)) +
  ylim(0, 4400) +
  geom_hline(yintercept = 3992, lty = "dashed", color = "darkgray") +
  annotate(geom = "text", y = 4000, x = 100, label = "3992 eQTLs with a single\ncredible causal variant")

### genotype x population interaction effect models ###

causal_sets_per_gene <- fread("~/Downloads/eQTL_finemapping.genes.txt.gz")

interaction <- fread("~/Downloads/fastQTL_eQTLs.maf-threshold.single-variant.superpop-interaction.txt.gz") %>%
  .[, gene_var := paste(geneID, variantID)] %>%
  .[, pval := F_pval]

interaction <- merge(interaction, causal_sets_per_gene, by = "geneID")

prop_by_threshold <- function(data, pval_threshold, n_causal_sets) {
  data.table(pval = pval_threshold,
             numCredibleSets = n_causal_sets,
             finemapped = nrow(data[numCredibleSets == n_causal_sets & pval < pval_threshold]),
             n_total_finemap = nrow(data[numCredibleSets == n_causal_sets])
  ) %>%
    return()
}

dt_interaction <- rbindlist(
  
  lapply(1:8,
         function(y) { rbindlist(
           
           
           pbmclapply(c(2:10 %o% 10^(-1:-7)), 
                      function(x) prop_by_threshold(interaction, x, y), 
                      mc.cores = getOption("mc.cores", 4L)
           )
         )
           
         })
  
)

dt_interaction[, group := "Single causal variant model"]

interaction_multivar <- fread("~/Downloads/susie_eQTLs.maf-threshold.multivariant.superpop-interaction.txt.gz") %>%
  .[, gene_var := paste(geneID, variantID)] %>%
  .[, pval := F_pval]

interaction_multivar <- merge(interaction_multivar, causal_sets_per_gene, by = "geneID")

prop_by_threshold <- function(data, pval_threshold, n_causal_sets) {
  data.table(pval = pval_threshold,
             numCredibleSets = n_causal_sets,
             finemapped = nrow(data[numCredibleSets == n_causal_sets & pval < pval_threshold]),
             n_total_finemap = nrow(data[numCredibleSets == n_causal_sets])
  ) %>%
    return()
}

dt_interaction_multivar <- rbindlist(
  
  lapply(1:8,
         function(y) { rbindlist(
           
           
           pbmclapply(c(2:10 %o% 10^(-1:-7)), 
                      function(x) prop_by_threshold(interaction_multivar, x, y), 
                      mc.cores = getOption("mc.cores", 4L)
           )
         )
           
         })
  
)

dt_interaction_multivar[, group := "Multiple causal variants model"]

dt <- rbind(dt_interaction, dt_interaction_multivar)

dt$group <- factor(dt$group, levels = c("Single causal variant model", "Multiple causal variants model"))

my_palette <- c('#ccebc5','#a8ddb5','#7bccc4','#43a2ca','#0868ac')

panel_interaction_left <- ggplot(data = dt[numCredibleSets %in% 1:5 & pval > 1e-6 &
                                             group == "Single causal variant model"], 
                                 aes(x = pval, y = finemapped / n_total_finemap, 
                                     color = factor(numCredibleSets, levels = 5:1))) +
  theme_bw() +
  geom_line() +
  geom_vline(xintercept = 5.969436e-06, lty = "dashed", color = "darkgray") +
  #scale_color_brewer(palette = "Set2", name = "") +
  xlab("p-value threshold") +
  ylab("Prop. significant interactions") +
  scale_x_log10() +
  facet_grid(. ~ group) +
  theme(legend.position = "none") +
  scale_color_manual(values = my_palette, name = "Number of\ncredible sets")

panel_interaction_right <- ggplot(data = dt[numCredibleSets %in% 1:5 & pval > 1e-6 &
                                              group == "Multiple causal variants model"], 
                                  aes(x = pval, y = finemapped / n_total_finemap, 
                                      color = factor(numCredibleSets, levels = 5:1))) +
  theme_bw() +
  geom_line() +
  geom_vline(xintercept = 4.052521e-06, lty = "dashed", color = "darkgray") +
  #scale_color_brewer(palette = "Set2", name = "") +
  xlab("p-value threshold") +
  ylab("Prop. significant interactions") +
  scale_x_log10() +
  facet_grid(. ~ group) +
  theme(legend.position = c(0.4, 0.57), legend.box.background = element_rect(colour = "darkgray")) +
  scale_color_manual(values = my_palette, name = "Number of credible sets")

plot_grid(panel_interaction_left, panel_interaction_right)

### stratify number of causal sets on pLI ###

dt <- fread("~/Downloads/eQTL_finemapping.genes.txt.gz") %>%
  setnames(., "geneID", "gene") %>%
  .[, gene := gsub("\\..*", "", gene)]

disp <- fread("~/Downloads/dispersions.csv")[expression_threshold.autosome == TRUE] %>%
  .[, gene := gsub("\\..*", "", gene)]

disp <- merge(dt, disp[expression_threshold.autosome == TRUE, c("gene", "baseMean", "dispersion")])

# restrict gnomad constraint data to MANE Select transcripts to avoid double-counting
plof <- fread("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/constraint/gnomad.v4.0.constraint_metrics.tsv") %>%
  setnames(., "gene", "gene_symbol") %>%
  .[mane_select == TRUE]

tx2gene <- rtracklayer::import("~/Downloads/gencode.v39.annotation.gtf.gz")
tx2gene <- data.table(tx2gene$transcript_id, tx2gene$gene_id) %>%
  setnames(., c("transcript", "gene")) %>%
  .[, gene := gsub("\\..*", "", gene)] %>%
  .[, transcript := gsub("\\..*", "", transcript)] %>%
  .[!is.na(transcript) & !is.na(gene) & !duplicated(transcript)]

plof <- merge(plof, tx2gene, by = "transcript")[, c("gene", "lof.pLI", "lof.oe_ci.upper")] %>%
  setnames(., "lof.pLI", "pLI")

disp <- merge(disp, plof, by = "gene", all.x = TRUE) %>%
  mutate(., pLI_decile = ntile(pLI, 10)) %>%
  as.data.table()

pli_dt <- group_by(disp[!is.na(pLI_decile)], numCredibleSets, pLI_decile == 10) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
pli_dt[, density := as.numeric(NA)]
pli_dt[constraint_criterion == TRUE, density := n / sum(n)]
pli_dt[constraint_criterion == FALSE, density := n / sum(n)]
pli_dt[, metric := "pLI"]

pli_p <- glm(data = disp[!is.na(pLI_decile)], 
             formula = numCredibleSets ~ baseMean + (pLI_decile == 10), 
             family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

pli_dt[, color_facet := as.character(NA)]
pli_dt[constraint_criterion == TRUE, color_facet := "Top 10%"]
pli_dt[constraint_criterion == FALSE, color_facet := "Lower 90%"]
pli_dt$color_facet <- factor(pli_dt$color_facet, levels = c("Top 10%", "Lower 90%"))

panel_cs_by_pli <- ggplot(data = pli_dt, aes(x = factor(numCredibleSets), y = density, fill = color_facet)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("Number of credible sets") +
  ylab("Proportion of eGenes") +
  scale_fill_manual(values = c("#CD7EAE","#668A99"), name = "pLI percentile") +
  theme(legend.position = c(0.7, 0.7), legend.box.background = element_rect(colour = "darkgray")) +
  annotate(y = 0.15, x = 8, label = paste("p =", formatC(pli_p, format = "e", digits = 2)), 
           color = "black", size = 4, geom = "text") +
  NULL

### stratify eQTL effect size by pLI ###

eqtl_effects <- fread("~/Downloads/kgpex.sample.filtered.covariate-corrected.ci95_wGT.aFC.txt.gz") %>%
  .[, gene_id := gsub("\\..*", "", gene_id)]

eqtl_effects[, pLI := gene_id %in% disp[pLI_decile == 10]$gene]
pli_eqtl_p <- wilcox.test(abs(eqtl_effects[pLI == TRUE]$log2_aFC), 
                          abs(eqtl_effects[pLI == FALSE]$log2_aFC))$p.value

panel_eqtl_beta_by_pli <- ggplot(data = eqtl_effects, 
                                 aes(x = abs(log2_aFC), fill = pLI, color = pLI)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  scale_color_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  theme(legend.position = "none") +
  ylab("Density")  +
  xlab(bquote('Absolute eQTL effect size (|' ~log[2]~(aFC)~'|)')) +
  # xlab(expression("Absolute eQTL effect size (|" ~ log_2_(aFC) ~ "|)")) +
  annotate(y = 1.5, x = 4, label = paste("p =", formatC(pli_eqtl_p, format = "e", digits = 2)), 
           color = "black", size = 4, geom = "text")

# test whether this holds for expression-increasing versus expression-decreasing alleles

ggplot(data = eqtl_effects, 
       aes(x = log2_aFC, fill = pLI, color = pLI)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  scale_color_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  theme(legend.position = "none") +
  ylab("Density")  +
  xlab(bquote('Absolute eQTL effect size (|' ~log[2]~(aFC)~'|)')) +
  # xlab(expression("Absolute eQTL effect size (|" ~ log_2_(aFC) ~ "|)")) +
  annotate(y = 1.5, x = 4, label = paste("p =", formatC(pli_eqtl_p, format = "e", digits = 2)), 
           color = "black", size = 4, geom = "text")

mean(abs(eqtl_effects[log2_aFC >0 & pLI == TRUE]$log2_aFC)) - mean(abs(eqtl_effects[log2_aFC >0 & pLI == FALSE]$log2_aFC))
wilcox.test(abs(eqtl_effects[log2_aFC >0 & pLI == TRUE]$log2_aFC), 
            abs(eqtl_effects[log2_aFC >0 & pLI == FALSE]$log2_aFC))$p.value

mean(abs(eqtl_effects[log2_aFC < 0 & pLI == TRUE]$log2_aFC)) - mean(abs(eqtl_effects[log2_aFC < 0 & pLI == FALSE]$log2_aFC))
wilcox.test(abs(eqtl_effects[log2_aFC < 0 & pLI == TRUE]$log2_aFC), 
            abs(eqtl_effects[log2_aFC < 0 & pLI == FALSE]$log2_aFC))$p.value

### plot multipanel ###

plot_grid(panel_cs_per_gene, panel_var_per_cs_log, 
          panel_cs_by_pli, panel_eqtl_beta_by_pli, 
          nrow = 2, align = "vh", axis = "lrtb",
          labels = paste0(LETTERS[1:6], "."))




