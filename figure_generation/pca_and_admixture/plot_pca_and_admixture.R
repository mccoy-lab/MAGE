library(data.table)
library(readxl)
library(ggh4x)
library(tidyverse)
library(BEDMatrix)
library(cowplot)
library(patchwork)
library(khroma)

cv_error <- fread("~/Downloads/cv_error.txt")

ggplot(data = cv_error, aes(x = V1, y = V2)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  xlab("k") +
  ylab("Cross-validation error") +
  xlim(0, 11)

meta <- fread("~/Downloads/igsr_samples.tsv") %>%
  .[, c("Sample name", "Population code", "Superpopulation code")] %>%
  setnames(., c("sid", "pop", "superpop")) %>%
  .[pop != ""]

kgpex_list <- fread("~/Downloads/selected_sample_map.txt", header = FALSE) %>%
  setnames(., c("proj_id", "sid"))

geuvadis_list <- fread("~/Downloads/igsr-geuvadis.tsv.tsv")[, 1] %>%
  setnames(., "sid")

afgr_list <- read_xlsx("~/Downloads/AFGR Supplementary Tables.xlsx", sheet = 1, skip = 1) %>%
  as.data.table() %>%
  setnames(., "1000G ID", "sid")

afgr_list2 <- fread("~/Downloads/AFGR.fam", header = FALSE) %>%
  setnames(., "V1", "sid")

gtex_list <- fread("https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") %>%
  .[SMAFRZE == "WGS"] %>%
  .[, sid := sub('^([^-]+-[^-]+).*', '\\1', SAMPID)]

eig <- fread("~/Downloads/plink.eigenvec") %>%
  setnames(., "V1", "sid") %>%
  merge(., meta, by = "sid", all.x = TRUE)

eig[grepl("GTEX", sid), sid := sub('^([^-]+-[^-]+).*', '\\1', sid)]

eig[sid %in% afgr_list2$sid, superpop := "AFR"]
eig[sid %in% afgr_list2$sid, pop := "MKK"]

ggplot(data = eig, aes(x = V3, y = V4, color = superpop))  +
  geom_point() +
  xlab("PC1") + 
  ylab("PC2")

ggplot(data = eig, aes(x = V3, y = V5, color = superpop))  +
  geom_point() +
  xlab("PC1") + 
  ylab("PC3")

ggplot(data = eig[superpop == "AFR"], aes(x = V3, y = V5, color = pop))  +
  geom_point() +
  xlab("PC1") + 
  ylab("PC3")

eig[V5 > 0.05]

###

dt_a <- eig[sid %in% gtex_list$sid] %>%
  .[, dataset := "GTEx"]

dt_b <- eig[sid %in% kgpex_list$sid]  %>%
  .[, dataset := "KGPEx"]

dt_c <- eig[sid %in% geuvadis_list$sid]  %>%
  .[, dataset := "Geuvadis"]

dt_d <- eig[sid %in% unique(c(afgr_list$sid, afgr_list2$sid))]  %>%
  .[, dataset := "AFGR"]

dt <- rbind(dt_a, dt_b, dt_c, dt_d)

dt$dataset <- factor(dt$dataset, levels = c("Geuvadis", "GTEx", "AFGR", "KGPEx"))

dt$superpop <- factor(dt$superpop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "NA"))

left_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "KGPEx"], aes(x = V3, y = V4), color = 'gray95', size = 0.5) +
  geom_point(data = dt[dataset == "KGPEx"], aes(x = V3, y = V4), color = '#e78ac3', size = 0.5) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "KGPEx", color = "black") +
  xlab("PC1") + ylab("PC2")

left_mid_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "Geuvadis"], aes(x = V3, y = V4), color = 'gray95', size = 0.5) +
  geom_point(data = dt[dataset == "Geuvadis"], aes(x = V3, y = V4), color = '#66c2a5', size = 0.5) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "Geuvadis", color = "black") +
  xlab("PC1") + ylab("PC2")

right_mid_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "GTEx"], aes(x = V3, y = V4), color = 'gray95', size = 0.5) +
  geom_point(data = dt[dataset == "GTEx"], aes(x = V3, y = V4), color = '#8da0cb', size = 0.5) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "GTEx", color = "black") +
  xlab("PC1") + ylab("PC2")

right_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "AFGR"], aes(x = V3, y = V4), color = 'gray95', size = 0.5) +
  geom_point(data = dt[dataset == "AFGR"], aes(x = V3, y = V4), color = '#fc8d62', size = 0.5) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "AFGR", color = "black") +
  xlab("PC1") + ylab("PC2")

pve <- fread("~/Downloads/pve.pca.eigenval") %>%
  .[, PC := 1:20] %>%
  .[, pve := (V1 / 2635.32) * 100]

#pve$PC <- factor(pve$PC, levels = paste0("PC", 1:20))

pve_panel <- ggplot(data = pve[PC < 11], aes(x = factor(PC), y = pve)) +
  geom_bar(stat = "identity") +
  xlab("Principal component") +
  ylab("Variance explained (%)") +
  theme_bw() +
  theme(panel.grid = element_blank())

#right_panel <- plot_grid(top_right_panel, middle_right_panel, bottom_right_panel, nrow = 3)

plot_grid(left_panel, left_mid_panel, right_mid_panel, right_panel, pve_panel, nrow = 1,
          rel_widths = c(1, 1, 1, 1, 0.9))



###

left_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "KGPEx"], aes(x = V3, y = V5), color = 'gray95', size = 1) +
  geom_point(data = dt[dataset == "KGPEx"], aes(x = V3, y = V5), color = '#e78ac3', size = 1) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "KGPEx", color = "black") +
  xlab("PC1") + ylab("PC3")

top_right_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "Geuvadis"], aes(x = V3, y = V5), color = 'gray95', size = 0.5) +
  geom_point(data = dt[dataset == "Geuvadis"], aes(x = V3, y = V5), color = '#66c2a5', size = 0.5) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "Geuvadis", color = "black") +
  xlab("PC1") + ylab("PC3")

middle_right_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "GTEx"], aes(x = V3, y = V5), color = 'gray95', size = 0.5) +
  geom_point(data = dt[dataset == "GTEx"], aes(x = V3, y = V5), color = '#8da0cb', size = 0.5) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "GTEx", color = "black") +
  xlab("PC1") + ylab("PC3")

bottom_right_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "AFGR"], aes(x = V3, y = V5), color = 'gray95', size = 0.5) +
  geom_point(data = dt[dataset == "AFGR"], aes(x = V3, y = V5), color = '#fc8d62', size = 0.5) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "AFGR", color = "black") +
  xlab("PC1") + ylab("PC3")

right_panel <- plot_grid(top_right_panel, middle_right_panel, bottom_right_panel, nrow = 3)

plot_grid(left_panel, right_panel, rel_widths = c(1, 0.35))

###

qmat <- fread("~/Downloads/merged_datasets.7.Q") %>%
  .[, sid := eig$sid]

qmat[grepl("GTEX", sid), sid := sub('^([^-]+-[^-]+).*', '\\1', sid)]

qmat <- qmat %>%
  .[, max_col := unname(unlist(lapply(1:nrow(qmat), function(x) which.max(qmat[x, 1:6]))))]

### try ordering on max ancestry

qmat_ordered <- rbind(
  qmat[max_col == 1] %>% setorder(., V1),
  qmat[max_col == 2] %>% setorder(., V2),
  qmat[max_col == 3] %>% setorder(., V3),
  qmat[max_col == 6] %>% setorder(., V6),
  qmat[max_col == 4] %>% setorder(., V4),
  qmat[max_col == 5] %>% setorder(., V5),
  qmat[max_col == 7] %>% setorder(., V7)
)

sid_ordered <- qmat_ordered$sid

my_palette <- c(
  "#77AADD", #V1
  "#BBCC33", #V2
  "#FFAABB", #V3
  "#44BB99", #V4
  "#EEDD88", #V5
  "#99DDFF", #V6
  "#EE8866"  #V7
)

qmat_melted <- melt(qmat[, 1:8], id.vars = "sid") %>%
  as.data.table() %>%
  merge(., meta, by = "sid", all.x = TRUE)

qmat_melted[sid %in% afgr_list2$sid, superpop := "AFR"]
qmat_melted[sid %in% afgr_list2$sid, pop := "MKK"]

qmat_melted$sid <- factor(qmat_melted$sid, levels = sid_ordered)

light <- colour("light")

panel_kgpex <- ggplot(data = qmat_melted[sid %in% kgpex_list$sid], aes(x = sid, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 1),
        panel.spacing = unit(0, "lines")) +
  facet_nested(. ~ superpop + pop, space = "free", scale = "free") +
  scale_fill_manual(values = my_palette) +
  xlab("") + ylab("Ancestry Prop.") +
  ggtitle("KGPEx (n = 731)")

panel_gtex <- ggplot(data = qmat_melted[sid %in% gtex_list$sid], aes(x = sid, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 1),
        panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values = my_palette) +
  xlab("") + ylab("Ancestry Prop.") +
  ggtitle("GTEx (n = 838; 4-706 RNA-seq samples with genotypes per tissue)")

panel_geuvadis <- ggplot(data = qmat_melted[sid %in% geuvadis_list$sid], aes(x = sid, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 1),
        panel.spacing = unit(0, "lines")) +
  facet_nested(. ~ superpop + pop, space = "free", scale = "free") +
  scale_fill_manual(values = my_palette) +
  xlab("") + ylab("Ancestry Prop.") +
  ggtitle("Geuvadis (n = 462)")

panel_afgr <- ggplot(data = qmat_melted[sid %in% unique(c(afgr_list$sid, afgr_list2$sid))], 
                     aes(x = sid, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 1),
        panel.spacing = unit(0, "lines")) +
  facet_nested(. ~ superpop + pop, space = "free", scale = "free") +
  scale_fill_manual(values = my_palette) +
  xlab("") + ylab("Ancestry Prop.") +
  ggtitle("AFGR (n = 599)")

plot_grid(panel_kgpex, NULL, panel_gtex, NULL, panel_geuvadis, NULL, panel_afgr,
          nrow = 7,
          rel_heights = c(1, -0.1, 0.8, -0.1, 1, -0.1, 1))



