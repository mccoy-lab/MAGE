library(data.table)
library(ggh4x)
library(tidyverse)
library(BEDMatrix)
library(cowplot)
library(patchwork)
library(khroma)

meta <- fread("~/Downloads/igsr_samples.tsv") %>%
  .[, c("Sample name", "Population code", "Superpopulation code")] %>%
  setnames(., c("sid", "pop", "superpop")) %>%
  .[pop != ""]

kgpex_list <- fread("~/Downloads/selected_sample_map.txt", header = FALSE) %>%
  setnames(., c("proj_id", "sid"))

geuvadis_list <- fread("~/Downloads/igsr-geuvadis.tsv.tsv")[, 1] %>%
  setnames(., "sid")

eig <- fread("~/Downloads/plink.eigenvec")

gtex_list <- data.table(sid = eig[grepl("GTEX", V1)]$V1)

###

dt_a <- eig[grepl("GTEX", V1)] %>%
  .[, dataset := "GTEx"]

dt_b <- eig[V1 %in% kgpex_list$sid]  %>%
  .[, dataset := "KGPEx"]

dt_c <- eig[V1 %in% geuvadis_list$sid]  %>%
  .[, dataset := "Geuvadis"]

dt <- rbind(dt_a, dt_b, dt_c)

dt$dataset <- factor(dt$dataset, levels = c("Geuvadis", "GTEx", "KGPEx"))

dt$superpop <- factor(dt$superpop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "NA"))

left_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "KGPEx"], aes(x = V3, y = V4), color = 'gray60', size = 1) +
  geom_point(data = dt[dataset == "KGPEx"], aes(x = V3, y = V4), color = '#e78ac3', size = 1) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "KGPEx", color = "black") +
  xlab("PC1") + ylab("PC2")

top_right_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "Geuvadis"], aes(x = V3, y = V4), color = 'gray60', size = 1) +
  geom_point(data = dt[dataset == "Geuvadis"], aes(x = V3, y = V4), color = '#66c2a5', size = 1) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "Geuvadis", color = "black") +
  xlab("PC1") + ylab("PC2")

bottom_right_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "GTEx"], aes(x = V3, y = V4), color = 'gray60', size = 1) +
  geom_point(data = dt[dataset == "GTEx"], aes(x = V3, y = V4), color = '#8da0cb', size = 1) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "GTEx", color = "black") +
  xlab("PC1") + ylab("PC2")

right_panel <- plot_grid(top_right_panel, bottom_right_panel, nrow = 2)

plot_grid(left_panel, right_panel, rel_widths = c(1, 0.55))

###

left_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "KGPEx"], aes(x = V3, y = V5), color = 'gray60', size = 1) +
  geom_point(data = dt[dataset == "KGPEx"], aes(x = V3, y = V5), color = '#e78ac3', size = 1) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "KGPEx", color = "black") +
  xlab("PC1") + ylab("PC3")

top_right_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "Geuvadis"], aes(x = V3, y = V5), color = 'gray60', size = 1) +
  geom_point(data = dt[dataset == "Geuvadis"], aes(x = V3, y = V5), color = '#66c2a5', size = 1) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "Geuvadis", color = "black") +
  xlab("PC1") + ylab("PC3")

bottom_right_panel <- ggplot() +
  theme_bw() +
  geom_point(data = dt[dataset != "GTEx"], aes(x = V3, y = V5), color = 'gray60', size = 1) +
  geom_point(data = dt[dataset == "GTEx"], aes(x = V3, y = V5), color = '#8da0cb', size = 1) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  annotate(geom = "text", x = 0.015, y = 0.03, size = 4, label = "GTEx", color = "black") +
  xlab("PC1") + ylab("PC3")

right_panel <- plot_grid(top_right_panel, bottom_right_panel, nrow = 2)

plot_grid(left_panel, right_panel, rel_widths = c(1, 0.55))

###

qmat <- fread("~/Downloads/1kgp_gtex_merged_maf5percent_LDindep.8.Q") %>%
  .[, sid := eig$V1] 

qmat <- qmat %>%
  .[, max_col := unname(unlist(lapply(1:nrow(qmat), function(x) which.max(qmat[x, 1:8]))))]

### try ordering on max ancestry

qmat_ordered <- rbind(
  qmat[max_col == 1] %>% setorder(., V1),
  qmat[max_col == 2] %>% setorder(., V2),
  qmat[max_col == 3] %>% setorder(., V3),
  qmat[max_col == 4] %>% setorder(., V4),
  qmat[max_col == 5] %>% setorder(., V5),
  qmat[max_col == 6] %>% setorder(., V6),
  qmat[max_col == 7] %>% setorder(., V7),
  qmat[max_col == 8] %>% setorder(., V8)
)

sid_ordered <- qmat_ordered$sid

qmat_melted <- melt(qmat[, 1:9], id.vars = "sid") %>%
  as.data.table() %>%
  merge(., meta, by = "sid", all.x = TRUE)

qmat_melted$sid <- factor(qmat_melted$sid, levels = sid_ordered)

light <- colour("light")

panel_kgpex <- ggplot(data = qmat_melted[sid %in% kgpex_list$sid], aes(x = sid, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 1),
        panel.spacing = unit(0, "lines")) +
  facet_nested(. ~ superpop + pop, space = "free", scale = "free") +
  scale_fill_manual(values = unname(light(8))) +
  xlab("") + ylab("Ancestry Prop.") +
  ggtitle("KGPEx (n = 731)")

panel_gtex <- ggplot(data = qmat_melted[sid %in% gtex_list$sid], aes(x = sid, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 1),
        panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values = unname(light(8))) +
  xlab("") + ylab("Ancestry Prop.") +
  ggtitle("GTEx (n = 953)")

panel_geuvadis <- ggplot(data = qmat_melted[sid %in% geuvadis_list$sid], aes(x = sid, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 1),
        panel.spacing = unit(0, "lines")) +
  facet_nested(. ~ superpop + pop, space = "free", scale = "free") +
  scale_fill_manual(values = unname(light(8))) +
  xlab("") + ylab("Ancestry Prop.") +
  ggtitle("Geuvadis (n = 465)")

plot_grid(panel_kgpex, panel_gtex, panel_geuvadis, 
          nrow = 3,
          rel_heights = c(1, 0.8, 1))

