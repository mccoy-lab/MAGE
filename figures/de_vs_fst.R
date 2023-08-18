library(data.table)
library(tidyverse)
library(cowplot)
library(broom)
library(khroma)

dt <- fread("~/Downloads/Fst_vs_DE.08042023.csv")

my_palette <- c(
  "#44BB99", 
  "#EEDD88",
  "#77AADD", 
  "#FFAABB",
  "#EE8866"
)

ggplot(data = dt[!is.na(PopSpecific_PadjDecile)], 
       aes(x = factor(PopSpecific_PadjDecile), y = mean_wcFst, fill = targetpop)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(. ~ targetpop) +
  xlab("Differential expression p-value decile") +
  ylab("Mean Fst of lead credible cis-eQTLs")
  
