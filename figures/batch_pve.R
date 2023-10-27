library(data.table)
library(tidyverse)
library(cowplot)
library(broom)
library(khroma)

dt_ex <- fread("~/Dropbox/projects/rna_seq_1kg/Replicate_ANOVA.csv") %>%
  .[, total_SS := SSRbatch + SSRsample + SSresiduals] %>%
  .[, Batch := SSRbatch/total_SS] %>%
  .[, Sample := SSRsample/total_SS] %>%
  melt(., id.vars = "Gene", measure.vars = c("Batch", "Sample")) %>%
  as.data.table()

panel_hist_ex <- ggplot(data = dt_ex, aes(x = value, fill = variable)) +
  geom_histogram(position = "identity", alpha = 0.6, binwidth = 0.0025) +
  xlab("Proportion variance explained") +
  ylab("Number of genes") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.7, 0.8), 
        legend.box.background = element_rect(colour = "darkgray"),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#CD7EAE","#668A99"), name = "")

dt_sp <- fread("~/Dropbox/projects/rna_seq_1kg/splicing_replicate_PVE.txt") %>%
  .[, Batch := batch_SS/total_SS] %>%
  .[, Sample := sample_SS/total_SS] %>%
  melt(., id.vars = "clusterID", measure.vars = c("Batch", "Sample")) %>%
  as.data.table()

panel_hist_sp <- ggplot(data = dt_sp, aes(x = value, fill = variable)) +
  geom_histogram(position = "identity", alpha = 0.6, binwidth = 0.0025) +
  xlab("Proportion variance explained") +
  ylab("Number of intron clusters") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.7, 0.8), 
        legend.box.background = element_rect(colour = "darkgray"),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#CD7EAE","#668A99"), name = "")

plot_grid(panel_hist_ex, panel_hist_sp,
          labels = c("A.", "B."))

