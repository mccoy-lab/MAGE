library(data.table)
library(tidyverse)
library(readxl)
library(broom)
library(cowplot)
library(GenomicRanges)
library(rtracklayer)

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

plof <- merge(plof, tx2gene, by = "transcript")[, c("gene", "lof.pLI", "lof.oe_ci.upper")]

disp <- merge(disp, plof, by = "gene", all.x = TRUE)

ds <- fread("~/Downloads/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz") %>%
  setnames(., c("gene_symbol", "pHaplo", "pTriplo"))

symbol2name <- fread("~/Downloads/gene_names.txt") %>%
  setnames(., c("gene", "gene_full", "gene_symbol"))

ds <- merge(ds, symbol2name, by = "gene_symbol")[, c("gene", "pHaplo", "pTriplo")]

disp <- merge(disp, ds, by = "gene", all.x = TRUE)

tx2gene <- fread("~/Downloads/mart_export.txt")[, c(1, 3)] %>%
  setnames(., c("gene", "Ensembl_transcript_id"))

hs <- fread("https://raw.githubusercontent.com/agarwal-i/loss-of-function-fitness-effects/main/out/Supplementary%20File%202.txt") %>%
  setorder(., -log10_map) %>%
  setnames(., "Gene", "gene_symbol") %>%
  merge(tx2gene, by = "Ensembl_transcript_id") %>%
  setnames(., "log10_map", "hs") %>%
  .[, hs := 10^hs]

disp <- merge(disp, hs[, c("gene", "hs")], by = "gene", all.x = TRUE)

rvis <- read_xlsx("~/Downloads/pgen.1003709.s002.xlsx") %>%
  as.data.table() %>%
  setnames(., c("gene_symbol", "rvis", "rvis_percentile"))
rvis <- merge(rvis, symbol2name, by = "gene_symbol")[, c("gene", "rvis")]

disp <- merge(disp, rvis, by = "gene", all.x = TRUE)

phylop_1000 <- fread("~/Downloads/kgpex.sample.filtered.inverse_normalized_tmm.all.tss.meanPhyloP.bed") %>%
  .[, 4:5] %>%
  setnames(., c("gene", "phylop_1000")) %>%
  .[, gene := gsub("\\..*", "", gene)]

disp <- merge(disp, phylop_1000, by = "gene", all.x = TRUE)

phylop_500 <- fread("~/Downloads/kgpex.sample.filtered.inverse_normalized_tmm.all.tss.broad_proximal_promoter.meanPhyloP.bed") %>%
  .[, 4:5] %>%
  setnames(., c("gene", "phylop_500")) %>%
  .[, gene := gsub("\\..*", "", gene)]

disp <- merge(disp, phylop_500, by = "gene", all.x = TRUE)

phylop_50 <- fread("~/Downloads/kgpex.sample.filtered.inverse_normalized_tmm.all.tss.proximal_promoter.meanPhyloP.bed") %>%
  .[, 4:5] %>%
  setnames(., c("gene", "phylop_50")) %>%
  .[, gene := gsub("\\..*", "", gene)]

disp <- merge(disp, phylop_50, by = "gene", all.x = TRUE)

#### restrict analysis to protein-coding genes

tx2gene <- rtracklayer::import("~/Downloads/gencode.v39.annotation.gtf.gz")
tx2gene <- data.table(tx2gene$gene_id, tx2gene$gene_type) %>%
  setnames(., c("gene", "gene_type")) %>%
  .[, gene := gsub("\\..*", "", gene)]

protein_coding <- tx2gene[gene_type == "protein_coding"]$gene
disp <- disp[gene %in% protein_coding]

# compute deciles of each metric
disp <- disp[!duplicated(gene)] %>%
  mutate(., expression_decile = ntile(baseMean, 10)) %>%
  mutate(., hs_decile = ntile(hs, 10)) %>%
  mutate(., rvis_decile = ntile(rvis, 10)) %>%
  mutate(., loeuf_decile = ntile(lof.oe_ci.upper, 10)) %>%
  mutate(., pLI_decile = ntile(lof.pLI, 10)) %>%
  mutate(., pHaplo_decile = ntile(pHaplo, 10)) %>%
  mutate(., pTriplo_decile = ntile(pTriplo, 10)) %>%
  mutate(., phylop_1000_decile = ntile(phylop_1000, 10)) %>%
  mutate(., phylop_500_decile = ntile(phylop_500, 10)) %>%
  mutate(., phylop_50_decile = ntile(phylop_50, 10)) %>%
  as.data.table()

# for each metric, define constrained set and compute densities in each credible set bin

loeuf_dt <- group_by(disp[!is.na(loeuf_decile)], numCredibleSets, loeuf_decile == 1) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
loeuf_dt[, density := as.numeric(NA)]
loeuf_dt[constraint_criterion == TRUE, density := n / sum(n)]
loeuf_dt[constraint_criterion == FALSE, density := n / sum(n)]
loeuf_dt[, metric := "LOEUF"]

# for each metric, compute a p-value comparing constrained and non-constrained sets

loeuf_p <- glm(data = disp[!is.na(loeuf_decile)], formula = numCredibleSets ~ baseMean + (loeuf_decile == 1), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

# repeat for other metrics

pli_dt <- group_by(disp[!is.na(pLI_decile)], numCredibleSets, pLI_decile == 10) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
pli_dt[, density := as.numeric(NA)]
pli_dt[constraint_criterion == TRUE, density := n / sum(n)]
pli_dt[constraint_criterion == FALSE, density := n / sum(n)]
pli_dt[, metric := "pLI"]

pli_p <- glm(data = disp[!is.na(pLI_decile)], formula = numCredibleSets ~ baseMean + (pLI_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

pHaplo_dt <- group_by(disp[!is.na(pHaplo_decile)], numCredibleSets, pHaplo_decile == 10) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
pHaplo_dt[, density := as.numeric(NA)]
pHaplo_dt[constraint_criterion == TRUE, density := n / sum(n)]
pHaplo_dt[constraint_criterion == FALSE, density := n / sum(n)]
pHaplo_dt[, metric := "pHaplo"]

pHaplo_p <- glm(data = disp[!is.na(pHaplo_decile)], formula = numCredibleSets ~ baseMean + (pHaplo_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

pTriplo_dt <- group_by(disp[!is.na(pTriplo_decile)], numCredibleSets, pTriplo_decile == 10) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
pTriplo_dt[, density := as.numeric(NA)]
pTriplo_dt[constraint_criterion == TRUE, density := n / sum(n)]
pTriplo_dt[constraint_criterion == FALSE, density := n / sum(n)]
pTriplo_dt[, metric := "pTriplo"]

pTriplo_p <- glm(data = disp[!is.na(pTriplo_decile)], formula = numCredibleSets ~ baseMean + (pTriplo_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

hs_dt <- group_by(disp[!is.na(hs_decile)], numCredibleSets, hs_decile == 10) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
hs_dt[, density := as.numeric(NA)]
hs_dt[constraint_criterion == TRUE, density := n / sum(n)]
hs_dt[constraint_criterion == FALSE, density := n / sum(n)]
hs_dt[, metric := "hs"]

hs_p <- glm(data = disp[!is.na(hs_decile)], formula = numCredibleSets ~ baseMean + (hs_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

rvis_dt <- group_by(disp[!is.na(rvis_decile)], numCredibleSets, rvis_decile == 1) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
rvis_dt[, density := as.numeric(NA)]
rvis_dt[constraint_criterion == TRUE, density := n / sum(n)]
rvis_dt[constraint_criterion == FALSE, density := n / sum(n)]
rvis_dt[, metric := "RVIS"]

rvis_p <- glm(data = disp[!is.na(rvis_decile)], formula = numCredibleSets ~ baseMean + (rvis_decile == 1), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

phylop_1000_dt <- group_by(disp[!is.na(phylop_1000_decile)], numCredibleSets, phylop_1000_decile == 10) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
phylop_1000_dt[, density := as.numeric(NA)]
phylop_1000_dt[constraint_criterion == TRUE, density := n / sum(n)]
phylop_1000_dt[constraint_criterion == FALSE, density := n / sum(n)]
phylop_1000_dt[, metric := "PhyloP ±1000 bp"]

phylop_1000_p <- glm(data = disp[!is.na(phylop_1000_decile)], formula = numCredibleSets ~ baseMean + (phylop_1000_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

phylop_500_dt <- group_by(disp[!is.na(phylop_500_decile)], numCredibleSets, phylop_500_decile == 10) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
phylop_500_dt[, density := as.numeric(NA)]
phylop_500_dt[constraint_criterion == TRUE, density := n / sum(n)]
phylop_500_dt[constraint_criterion == FALSE, density := n / sum(n)]
phylop_500_dt[, metric := "PhyloP [-500, 0] bp"]

phylop_500_p <- glm(data = disp[!is.na(phylop_500_decile)], formula = numCredibleSets ~ baseMean + (phylop_500_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

phylop_50_dt <- group_by(disp[!is.na(phylop_50_decile)], numCredibleSets, phylop_50_decile == 10) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
phylop_50_dt[, density := as.numeric(NA)]
phylop_50_dt[constraint_criterion == TRUE, density := n / sum(n)]
phylop_50_dt[constraint_criterion == FALSE, density := n / sum(n)]
phylop_50_dt[, metric := "PhyloP [-50, 0] bp"]

phylop_50_p <- glm(data = disp[!is.na(phylop_50_decile)], formula = numCredibleSets ~ baseMean + (phylop_50_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

# combine all metrics for plotting

constraint_dt <- rbind(pli_dt, loeuf_dt, pHaplo_dt, pTriplo_dt, hs_dt, rvis_dt, phylop_1000_dt, phylop_500_dt, phylop_50_dt)
constraint_dt[, density_flip := density]
constraint_dt[constraint_criterion == FALSE, density_flip := -1 * density]
constraint_dt[, y_facet := as.character(NA)]
constraint_dt[constraint_criterion == TRUE, y_facet := "Top 10%"]
constraint_dt[constraint_criterion == FALSE, y_facet := "Lower 90%"]

constraint_dt$y_facet <- factor(constraint_dt$y_facet, levels = c("Top 10%", "Lower 90%"))

constraint_dt$metric <- factor(constraint_dt$metric, 
                               levels = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS", "PhyloP ±1000 bp", "PhyloP [-500, 0] bp", "PhyloP [-50, 0] bp"))

pval_dt <- data.table(metric = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS", "PhyloP ±1000 bp", "PhyloP [-500, 0] bp", "PhyloP [-50, 0] bp"),
                      pval = c(pli_p, loeuf_p, pHaplo_p, pTriplo_p, hs_p, rvis_p, phylop_1000_p, phylop_500_p, phylop_50_p),
                      y_facet = "Top 10%")

pval_dt$metric <- factor(pval_dt$metric, 
                         levels = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS", "PhyloP ±1000 bp", "PhyloP [-500, 0] bp", "PhyloP [-50, 0] bp"))

top_panel_gene <- ggplot(data = constraint_dt[!grepl("PhyloP", metric)], 
                         aes(x = numCredibleSets, y = density, fill = y_facet)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlim(-1, 11) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  facet_grid(. ~ metric, scales = "free", space = "free") + 
  ylab("Proportion of eGenes") +
  xlab("Number of credible sets") +
  scale_fill_manual(values = c("#CD7EAE","#668A99"), name = "") +
  geom_text(data = pval_dt[!grepl("PhyloP", metric)], aes(x = 7, y = 0.35, 
                                                          label = paste("p =", formatC(pval, format = "e", digits = 2))), 
            color = "black", size = 3)

top_panel_promoter <- ggplot(data = constraint_dt[grepl("PhyloP", metric)], 
                             aes(x = numCredibleSets, y = density, fill = y_facet)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlim(-1, 11) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  facet_grid(. ~ metric, scales = "free", space = "free") + 
  ylab("Proportion of eGenes") +
  xlab("Number of credible sets") +
  scale_fill_manual(values = c("#CD7EAE","#668A99"), name = "") +
  geom_text(data = pval_dt[grepl("PhyloP", metric)], aes(x = 7, y = 0.35, 
                                                         label = paste("p =", formatC(pval, format = "e", digits = 2))), 
            color = "black", size = 3)

#### contrast eQTL effect sizes between constrained and non-constrained sets

eqtl_effects <- fread("~/Downloads/kgpex.sample.filtered.covariate-corrected.ci95_wGT.aFC.txt.gz") %>%
  .[, gene_id := gsub("\\..*", "", gene_id)]

eqtl_effects[, pLI := gene_id %in% disp[pLI_decile == 10]$gene]
eqtl_effects[, LOEUF := gene_id %in% disp[loeuf_decile == 1]$gene]
eqtl_effects[, pHaplo := gene_id %in% disp[pHaplo_decile == 10]$gene]
eqtl_effects[, pTriplo := gene_id %in% disp[pTriplo_decile == 10]$gene]
eqtl_effects[, hs := gene_id %in% disp[hs_decile == 10]$gene]
eqtl_effects[, RVIS := gene_id %in% disp[rvis_decile == 1]$gene]
eqtl_effects[, `PhyloP ±1000 bp` := gene_id %in% disp[phylop_1000_decile == 10]$gene]
eqtl_effects[, `PhyloP [-500, 0] bp` := gene_id %in% disp[phylop_500_decile == 10]$gene]
eqtl_effects[, `PhyloP [-50, 0] bp` := gene_id %in% disp[phylop_50_decile == 10]$gene]

eqtl_effects_melted <- melt(eqtl_effects, id.vars = c("variant_id", "gene_id", "log2_aFC"), 
                            measure.vars = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS", "PhyloP ±1000 bp", "PhyloP [-500, 0] bp", "PhyloP [-50, 0] bp"))

eqtl_effects_melted[, is_constrained := as.character("NA")]
eqtl_effects_melted[value == TRUE, is_constrained := "Top 10%"]
eqtl_effects_melted[value == FALSE, is_constrained := "Lower 90%"]

pli_eqtl_p <- wilcox.test(log(abs(eqtl_effects[pLI == TRUE]$log2_aFC)), log(abs(eqtl_effects[pLI == FALSE]$log2_aFC)))$p.value
loeuf_eqtl_p <- wilcox.test(log(abs(eqtl_effects[LOEUF == TRUE]$log2_aFC)), log(abs(eqtl_effects[LOEUF == FALSE]$log2_aFC)))$p.value
pHaplo_eqtl_p <- wilcox.test(log(abs(eqtl_effects[pHaplo == TRUE]$log2_aFC)), log(abs(eqtl_effects[pHaplo == FALSE]$log2_aFC)))$p.value
pTriplo_eqtl_p <- wilcox.test(log(abs(eqtl_effects[pTriplo == TRUE]$log2_aFC)), log(abs(eqtl_effects[pTriplo == FALSE]$log2_aFC)))$p.value
hs_eqtl_p <- wilcox.test(log(abs(eqtl_effects[hs == TRUE]$log2_aFC)), log(abs(eqtl_effects[hs == FALSE]$log2_aFC)))$p.value
rvis_eqtl_p <- wilcox.test(log(abs(eqtl_effects[RVIS == TRUE]$log2_aFC)), log(abs(eqtl_effects[RVIS == FALSE]$log2_aFC)))$p.value
phylop_1000_eqtl_p <- wilcox.test(log(abs(eqtl_effects[`PhyloP ±1000 bp` == TRUE]$log2_aFC)), log(abs(eqtl_effects[`PhyloP ±1000 bp` == FALSE]$log2_aFC)))$p.value
phylop_500_eqtl_p <- wilcox.test(log(abs(eqtl_effects[`PhyloP [-500, 0] bp` == TRUE]$log2_aFC)), log(abs(eqtl_effects[`PhyloP [-500, 0] bp` == FALSE]$log2_aFC)))$p.value
phylop_50_eqtl_p <- wilcox.test(log(abs(eqtl_effects[`PhyloP [-50, 0] bp` == TRUE]$log2_aFC)), log(abs(eqtl_effects[`PhyloP [-50, 0] bp` == FALSE]$log2_aFC)))$p.value

pval_eqtl_dt <- data.table(variable = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS", "PhyloP ±1000 bp", "PhyloP [-500, 0] bp", "PhyloP [-50, 0] bp"),
                           pval = c(pli_eqtl_p, loeuf_eqtl_p, pHaplo_eqtl_p, pTriplo_eqtl_p, hs_eqtl_p, rvis_eqtl_p, phylop_1000_eqtl_p, phylop_500_eqtl_p, phylop_50_eqtl_p),
                           is_constrained = "Top 10%")

pval_eqtl_dt$variable <- factor(pval_eqtl_dt$variable, 
                                levels = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS", "PhyloP ±1000 bp", "PhyloP [-500, 0] bp", "PhyloP [-50, 0] bp"))

eqtl_effects_melted$variable <- factor(eqtl_effects_melted$variable, 
                                       levels = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS", "PhyloP ±1000 bp", "PhyloP [-500, 0] bp", "PhyloP [-50, 0] bp"))

lower_panel_gene <- ggplot(data = eqtl_effects_melted[!grepl("PhyloP", variable)], 
                           aes(x = abs(log2_aFC), fill = is_constrained)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  facet_grid(. ~ variable, scales = "free", space = "free") + 
  scale_fill_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  scale_color_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  theme(legend.position = "none") +
  ylab("Density")  +
  xlab(bquote('Absolute eQTL effect size (|' ~log[2]~(aFC)~'|)')) +
  geom_text(data = pval_eqtl_dt[!grepl("PhyloP", variable)], aes(x = 4, y = 2.8, 
                                                                 label = paste("p =", formatC(pval, format = "e", digits = 2))), 
            color = "black", size = 3)

lower_panel_promoter <- ggplot(data = eqtl_effects_melted[grepl("PhyloP", variable)], 
                               aes(x = abs(log2_aFC), fill = is_constrained)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  facet_grid(. ~ variable, scales = "free", space = "free") + 
  scale_fill_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  scale_color_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  theme(legend.position = "none") +
  ylab("Density")  +
  xlab(bquote('Absolute eQTL effect size (|' ~log[2]~(aFC)~'|)')) +
  geom_text(data = pval_eqtl_dt[grepl("PhyloP", variable)], aes(x = 4, y = 2.8, 
                                                                label = paste("p =", formatC(pval, format = "e", digits = 2))), 
            color = "black", size = 3)

### multipanel grid ###

gene_panel <- plot_grid(top_panel_gene + ggtitle("Gene-level metrics"), lower_panel_gene, 
                        nrow = 2, align = "v", axis = "lr", labels = c("A.", ""))

promoter_panel <- plot_grid(top_panel_promoter + ggtitle("Promoter-level metrics"), lower_panel_promoter,
                            nrow = 2, align = "v", axis = "lr", labels = c("B.", ""))

plot_grid(gene_panel, promoter_panel, rel_widths = c(1, 0.5), nrow = 2, align = "v", axis = "l")



