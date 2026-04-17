#!/usr/bin/env Rscript
# Instructions from:
# https://github.com/davidaknowles/leafcutter/issues/265#issuecomment-2679783289

install.packages("BiocManager",
                 dependencies = TRUE,
                 repos = "https://cloud.r-project.org")
BiocManager::install("Biobase",
                     ask = FALSE,
                     update = FALSE)
install.packages("TailRank",
                 dependencies = TRUE,
                 repos = "https://cloud.r-project.org")
BiocManager::install("DirichletMultinomial", 
                     ask = FALSE,
                     update = FALSE)
devtools::install_github("davidaknowles/leafcutter/leafcutter",
                         ref = "psi_2019",
                         dependencies = TRUE,
                         upgrade = "never")

# Print session info to a file
library(leafcutter)
sink("results/leafcutter_R_loaded_packages.txt")
print(sessionInfo())
sink()