library(data.table)
library(tidyverse)
library(cowplot)
library(broom)
library(khroma)

dt_ex <- fread("~/Downloads/PVE.expression.filtered_genes.csv") %>%
  .[, mean_superpop_plus_subpop := Mean_superpop_proportions + Mean_super.subpop_proportions] %>%
  .[, c("V1", "Mean_superpop_proportions", "mean_superpop_plus_subpop")] %>%
  melt(id.vars = "V1") 

dt_ex[variable == "Mean_superpop_proportions", variable := "Continental group"]
dt_ex[variable == "mean_superpop_plus_subpop", variable := "Population"]

panel_pve_expression <- ggplot(data = dt_ex, aes(x = value, fill = variable)) +
  geom_histogram(position = "identity", alpha = 0.6, binwidth = 0.0025) +
  xlab("Proportion variance explained") +
  ylab("Number of genes") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.6, 0.7), 
        legend.box.background = element_rect(colour = "darkgray"),
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#CD7EAE","#668A99"), name = "") +
  xlim(0, 0.5) +
  facet_grid(. ~ "Gene expression level")

dt_sp <- fread("~/Downloads/kgpex.sample.filtered.anova-ss.txt") %>%
  .[, superpop_PVE := superpop_SS / batch_sex_residual_SS] %>%
  .[, subpop_PVE := subpop_SS / batch_sex_residual_SS] %>%
  .[, c("clusterID", "superpop_PVE", "subpop_PVE")] %>%
  melt(id.vars = "clusterID") 

dt_sp[variable == "superpop_PVE", variable := "Continental group"]
dt_sp[variable == "subpop_PVE", variable := "Population"]

panel_pve_splicing <- ggplot(data = dt_sp, aes(x = value, fill = variable)) +
  geom_histogram(position = "identity", alpha = 0.6, binwidth = 0.0025) +
  xlab("Proportion variance explained") +
  ylab("Number of splicing clusters") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_manual(values = c("#CD7EAE","#668A99")) +
  xlim(0, 0.5) +
  facet_grid(. ~ "Alternative splicing")

### variance within population ###

my_palette <- c(
  "#44BB99", 
  "#FFAABB",
  "#EE8866",
  "#77AADD", 
  "#EEDD88"
)

dt_ex_resid <- fread("~/Downloads/PVE_ANOVA.Superpop_ResidualVariance.csv")

dt_ex_resid$subset <- factor(dt_ex_resid$subset, levels = c("AFR", "EUR", "SAS", "EAS", "AMR"))

var_by_superpop <- lm(data = dt_ex_resid, formula = log10(variance) ~ -1 + subset) %>%
  summary() %>%
  tidy() %>%
  mutate(term=str_replace(term,"subset",""))

var_by_superpop$term <- factor(var_by_superpop$term, 
                               levels = c("AFR", "EUR", "EAS", "SAS", "AMR"))

var_by_superpop <- group_by(dt_ex_resid, subset) %>%
  summarize(median_var = median(variance)) %>%
  as.data.table()

panel_resid_var_expression <- ggplot(data = dt_ex_resid, aes(x = subset, y = variance)) +
  geom_violin(aes(fill = subset), color = NA) +
  geom_boxplot(width = 0.25, outlier.alpha = 0) +
  scale_y_log10(limits = c(0.01, 10), labels = function(x) format(x, scientific = TRUE)) +
  theme_bw() +
  scale_fill_manual(values = my_palette) +
  theme(legend.position = "none") +
  annotation_logticks(sides = "l") +
  ylab("Variance") +
  xlab("Continental group") +
  facet_grid(. ~ "Gene expression level") 

dt_sp_resid <- fread("~/Downloads/kgpex.sample.filtered.superpop-variance.txt") %>%
  melt(id.vars = "clusterID") %>%
  as.data.table()

dt_sp_resid$variable <- factor(dt_sp_resid$variable, 
                               levels = c("AFR", "EUR", "EAS", "SAS", "AMR"))

# this is a placeholder - replace with actual splicing data 
panel_resid_var_splicing <- ggplot(data = dt_sp_resid, aes(x = variable, y = value)) +
  geom_violin(aes(fill = variable), color = NA) +
  geom_boxplot(width = 0.25, outlier.alpha = 0) +
  scale_y_log10(limits = c(1e-4, 1), labels = function(x) format(x, scientific = TRUE)) +
  theme_bw() +
  scale_fill_manual(values = my_palette) +
  theme(legend.position = "none") +
  annotation_logticks(sides = "l") +
  ylab("Variance") +
  xlab("Continental group") +
  facet_grid(. ~ "Alternative splicing") 
  
plot_grid(panel_pve_expression, panel_resid_var_expression, 
          panel_pve_splicing, panel_resid_var_splicing, 
          labels = c("A.", "B.", "C.", "D."), ncol = 2,
          align = "v", axis = "lr")
