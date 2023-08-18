library(data.table)
library(tidyverse)
library(readxl)
library(broom)
library(cowplot)

dt <- fread("~/Downloads/eQTL_finemapping.genes.txt.gz") %>%
  setnames(., "geneID", "gene") %>%
  .[, gene := gsub("\\..*", "", gene)]

disp <- fread("~/Downloads/dispersions.csv")[expression_threshold.autosome == TRUE] %>%
  .[, gene := gsub("\\..*", "", gene)]

disp <- merge(dt, disp[expression_threshold.autosome == TRUE, c("gene", "baseMean", "dispersion")])
  
plof <- fread("bgzip -c -d https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz") %>%
  setnames(., "gene", "gene_symbol") %>%
  setnames(., "gene_id", "gene")

disp <- merge(disp, plof, by = "gene", all.x = TRUE)

ds <- fread("~/Downloads/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz") %>%
  setnames(., c("gene_symbol", "pHaplo", "pTriplo"))

symbol2name <- fread("~/Downloads/gene_names.txt") %>%
  setnames(., c("gene", "gene_full", "gene_symbol"))

ds <- merge(ds, symbol2name, by = "gene_symbol")

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

disp <- disp[!duplicated(gene)] %>%
  mutate(., expression_decile = ntile(baseMean, 10)) %>%
  mutate(., hs_decile = ntile(hs, 10)) %>%
  mutate(., rvis_decile = ntile(rvis, 10)) %>%
  mutate(., loeuf_decile = ntile(oe_lof_upper, 10)) %>%
  mutate(., pLI_decile = ntile(pLI, 10)) %>%
  mutate(., pHaplo_decile = ntile(pHaplo, 10)) %>%
  mutate(., pTriplo_decile = ntile(pTriplo, 10)) %>%
  as.data.table()

loeuf_dt <- group_by(disp[!is.na(loeuf_decile)], numCredibleSets, loeuf_decile == 1) %>%
  summarize(., n = n()) %>%
  as.data.table() %>%
  setnames(., c("numCredibleSets", "constraint_criterion", "n")) %>%
  .[!(is.na(constraint_criterion))]
loeuf_dt[, density := as.numeric(NA)]
loeuf_dt[constraint_criterion == TRUE, density := n / sum(n)]
loeuf_dt[constraint_criterion == FALSE, density := n / sum(n)]
loeuf_dt[, metric := "LOEUF"]

loeuf_p <- glm(data = disp[!is.na(loeuf_decile)], formula = numCredibleSets ~ baseMean + (loeuf_decile == 1), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[3,] %>%
  .$p.value

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

pTriplo_p <- glm(data = disp[!is.na(pTriplo_decile)], formula = numCredibleSets ~ (pTriplo_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[2,] %>%
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

hs_p <- glm(data = disp[!is.na(hs_decile)], formula = numCredibleSets ~ (hs_decile == 10), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[2,] %>%
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

rvis_p <- glm(data = disp[!is.na(rvis_decile)], formula = numCredibleSets ~ (rvis_decile == 1), family = "quasipoisson") %>%
  tidy() %>%
  as.data.table() %>%
  .[2,] %>%
  .$p.value

constraint_dt <- rbind(pli_dt, loeuf_dt, pHaplo_dt, pTriplo_dt, hs_dt, rvis_dt)
constraint_dt[, density_flip := density]
constraint_dt[constraint_criterion == FALSE, density_flip := -1 * density]
constraint_dt[, y_facet := as.character(NA)]
constraint_dt[constraint_criterion == TRUE, y_facet := "Top 10%"]
constraint_dt[constraint_criterion == FALSE, y_facet := "Lower 90%"]

constraint_dt$y_facet <- factor(constraint_dt$y_facet, levels = c("Top 10%", "Lower 90%"))

constraint_dt$metric <- factor(constraint_dt$metric, 
                               levels = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS"))

pval_dt <- data.table(metric = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS"),
           pval = c(pli_p, loeuf_p, pHaplo_p, pTriplo_p, hs_p, rvis_p),
           y_facet = "Top 10%")

pval_dt$metric <- factor(pval_dt$metric, 
                         levels = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS"))

top_panel <- ggplot(data = constraint_dt, aes(x = numCredibleSets, y = density, fill = y_facet)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlim(-1, 11) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  facet_grid(. ~ metric, scales = "free") + 
  ylab("Proportion of eGenes") +
  xlab("Number of credible sets") +
  scale_fill_manual(values = c("#CD7EAE","#668A99"), name = "") +
  geom_text(data = pval_dt, aes(x = 7, y = 0.35, 
                                label = paste("p =", formatC(pval, format = "e", digits = 2))), 
            color = "black", size = 3)

####

eqtl_effects <- fread("~/Downloads/eQTL_finemapping.credibleSets.marginalEffects.txt.gz") %>%
  .[, gene_id := gsub("\\..*", "", geneID)]

eqtl_effects[, pLI := gene_id %in% disp[pLI_decile == 10]$gene]
eqtl_effects[, LOEUF := gene_id %in% disp[loeuf_decile == 1]$gene]
eqtl_effects[, pHaplo := gene_id %in% disp[pHaplo_decile == 10]$gene]
eqtl_effects[, pTriplo := gene_id %in% disp[pTriplo_decile == 10]$gene]
eqtl_effects[, hs := gene_id %in% disp[hs_decile == 10]$gene]
eqtl_effects[, RVIS := gene_id %in% disp[rvis_decile == 1]$gene]

eqtl_effects_melted <- melt(eqtl_effects, id.vars = c("topHitVariantID", "gene_id", "topHitMarginalSlope"), 
                        measure.vars = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS"))

eqtl_effects_melted[, is_constrained := as.character("NA")]
eqtl_effects_melted[value == TRUE, is_constrained := "Top 10%"]
eqtl_effects_melted[value == FALSE, is_constrained := "Lower 90%"]

pli_eqtl_p <- t.test(log(abs(eqtl_effects[pLI == TRUE]$topHitMarginalSlope)), log(abs(eqtl_effects[pLI == FALSE]$topHitMarginalSlope)))$p.value
loeuf_eqtl_p <- t.test(log(abs(eqtl_effects[LOEUF == TRUE]$topHitMarginalSlope)), log(abs(eqtl_effects[LOEUF == FALSE]$topHitMarginalSlope)))$p.value
pHaplo_eqtl_p <- t.test(log(abs(eqtl_effects[pHaplo == TRUE]$topHitMarginalSlope)), log(abs(eqtl_effects[pHaplo == FALSE]$topHitMarginalSlope)))$p.value
pTriplo_eqtl_p <- t.test(log(abs(eqtl_effects[pTriplo == TRUE]$topHitMarginalSlope)), log(abs(eqtl_effects[pTriplo == FALSE]$topHitMarginalSlope)))$p.value
hs_eqtl_p <- t.test(log(abs(eqtl_effects[hs == TRUE]$topHitMarginalSlope)), log(abs(eqtl_effects[hs == FALSE]$topHitMarginalSlope)))$p.value
rvis_eqtl_p <- t.test(log(abs(eqtl_effects[RVIS == TRUE]$topHitMarginalSlope)), log(abs(eqtl_effects[RVIS == FALSE]$topHitMarginalSlope)))$p.value

pval_eqtl_dt <- data.table(variable = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS"),
                      pval = c(pli_eqtl_p, loeuf_eqtl_p, pHaplo_eqtl_p, pTriplo_eqtl_p, hs_eqtl_p, rvis_eqtl_p),
                      is_constrained = "Top 10%")

eqtl_effects_melted$variable <- factor(eqtl_effects_melted$variable, 
                         levels = c("pLI", "LOEUF", "pHaplo", "pTriplo", "hs", "RVIS"))

lower_panel <- ggplot(data = eqtl_effects_melted, aes(x = abs(topHitMarginalSlope), fill = is_constrained)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  facet_grid(. ~ variable, scales = "free") + 
  scale_fill_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  scale_color_manual(values = rev(c("#CD7EAE","#668A99")), name = "") +
  theme(legend.position = "none") +
  ylab("Density")  +
  xlab(expression("Absolute eQTL effect size (|" ~ hat(beta) ~ "|)")) +
  geom_text(data = pval_eqtl_dt, aes(x = 1.7, y = 2.8, 
                                     label = paste("p =", formatC(pval, format = "e", digits = 2))), 
            color = "black", size = 3)

### multipanel grid ###

plot_grid(top_panel, lower_panel, nrow = 2, align = "v", axis = "lr")
