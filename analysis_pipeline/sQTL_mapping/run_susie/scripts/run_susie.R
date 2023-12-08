#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(Rfast)
  library(Matrix)
  library(vcfR)
  library(susieR)
  library(argparser)
})
set.seed(1)


#==================#
# Define functions #
#==================#

# This function regresses covariates out of both the response variable (y, expression)
# and explanatory variables (X, genotypes). It is taken from the SuSiE wiki:
# https://stephenslab.github.io/susieR/articles/finemapping.html#a-note-on-covariate-adjustment
remove.covariate.effects <- function (X, Z, y) {
  # include the intercept term
  if (any(Z[,1]!=1)) Z = cbind(1, Z)
  A   <- forceSymmetric(crossprod(Z))
  SZy <- as.vector(solve(A,c(y %*% Z)))
  SZX <- as.matrix(solve(A,t(Z) %*% X))
  y <- y - c(Z %*% SZy)
  X <- X - Z %*% SZX
  return(list(X=X, y=y, SZy=SZy, SZX=SZX))
}

numeric.matrix <- function(Z) {
  for (cov in colnames(Z)) {
    if (suppressWarnings(all(!is.na(as.numeric(Z[,cov]))))) {
      Z[,cov] <- as.numeric(Z[,cov])
    } else {
      Z[,cov] <- as.factor(Z[,cov])
    }
  }
  return(model.matrix(~., data=Z))
}


#===================#
# Collect arguments #
#===================#

p <- arg_parser("Run SuSiE sQTL fine-mapping for a single intron")
p <- add_argument(p, '--intronID', short='-g', help='ID of intron to run fine-mapping for')
p <- add_argument(p, '--gtPrefix', short='-m', help='Prefix (including path) of genotype matrix files. There should be three files, ending ".012", ".012.indv", and ".012.snps"')
p <- add_argument(p, '--usageBed', short='-e', help='Bed file with normalized intron usage values for each sample')
p <- add_argument(p, '--nominalResults', short='-n', help='Nominal results file from FastQTL')
p <- add_argument(p, '--covFile', short='-c', help='File with additional covariates to be regressed out')
p <- add_argument(p, '--outDir', short='-o', help="Directory to write output files to")
args <- parse_args(p)


#==============#
# Prepare Data #
#==============#

cat('\nRunning SuSiE for intron: ', args$intronID, '\n')
cat('\t- Reading in genotypes\n')
snp_ids <- read.table(paste0(args$gtPrefix, ".012.snps"))[,1]
sample_ids <- read.table(paste0(args$gtPrefix, ".012.indv"))[,1]
X <- as.matrix(read.table(paste0(args$gtPrefix, ".012"), sep='\t', header=FALSE, row.names=1))
rownames(X) <- sample_ids
colnames(X) <- snp_ids


cat('\t- Reading in expression data\n')
expression_data <- read.table(args$usageBed, sep='\t', header=TRUE, comment.char="")
expression_data <- expression_data[,4:ncol(expression_data)]
Y <- as.numeric(expression_data[expression_data[,1]==args$intronID, sample_ids])
names(Y) <- sample_ids


cat('\t- Reading in covariates\n')
Z <- as.data.frame(t(read.table(args$covFile, sep='\t', header=FALSE)))
colnames(Z) <- Z[1,]
rownames(Z) <- Z[,1]
Z <- Z[-1,-1]
Z <- Z[sample_ids,]
Z <- numeric.matrix(Z)


cat('\t- Regressing out covariates\n')
out = remove.covariate.effects(X, Z, Y)
X <- out$X
Y <- out$y


cat('\t- Reading in nominal results\n')
target_gene_results <- read.table(args$nominalResults, sep='\t', header=TRUE)
target_gene_results <- arrange(target_gene_results, factor(variant_id, levels=snp_ids))
z_scores <- target_gene_results$slope / target_gene_results$slope_se


cat('\t- Calculating variances\n')
gene_ld_matrix <- cor(X)
target_gene_var <- var(Y)
n <- length(Y)


#=================#
# Fit SuSiE Model #
#=================#

cat('\t- Fitting with Susie\n')
fitted_results <- susie_rss(
  z=z_scores,
  n=n,
  R=gene_ld_matrix,
  var_y=target_gene_var,  
  L=10,
  coverage=0.95,
  min_abs_corr = 0.5,
  estimate_residual_variance=TRUE
)

# fitted_results <- susie(
#   X,
#   Y,
#   L=10,
#   coverage=0.95,
#   min_abs_corr = 0.5,
#   estimate_residual_variance=TRUE
# )

#===================#
# Save Model Output #
#===================#

cat('\t- Writing output\n')
snp_output_file <- paste0(args$outDir, "/", args$intronID, "_full_pip_results.txt")
cs_output_file <- paste0(args$outDir, "/", args$intronID, "_cs_results.txt")

results <- data.frame(variant_id=character(), pip=numeric(), credible_set=character())

for (variant_idx in 1:length(snp_ids)) {
  results[nrow(results)+1,] = c(snp_ids[variant_idx], fitted_results$pip[variant_idx], NA)
}

cs_results <- data.frame(credible_set=character(), cs_nvars=integer(), cs_coverage=numeric(), cs_min_corr=numeric(), cs_mean_corr=numeric(), cs_median_corr=numeric())

if (length(names(fitted_results$set$cs)) > 0) {
  for (i in 1:length(names(fitted_results$set$cs))) {
    cred_set <- names(fitted_results$set$cs)[i]
    cs_coverage <- fitted_results$set$coverage[i]
    cs_nvars <- length(fitted_results$set$cs[[i]])
    min_corr <- fitted_results$set$purity[cred_set, 'min.abs.corr']
    mean_corr <- fitted_results$set$purity[cred_set, 'mean.abs.corr']
    median_corr <- fitted_results$set$purity[cred_set, 'median.abs.corr']
    cs_results[nrow(cs_results)+1,] = c(cred_set, cs_nvars, cs_coverage, min_corr, mean_corr, median_corr)
    for (var_index in fitted_results$set$cs[cred_set][[1]]) {
      results[var_index, 'credible_set'] = cred_set
    }
  }
}

write.table(results, snp_output_file, append=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(cs_results, cs_output_file, append=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

cat('Done.\n\n')
warnings()
