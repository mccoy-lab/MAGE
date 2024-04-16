#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(Rfast)
  library(Matrix)
  library(vcfR)
  library(susieR)
  library(coloc)
  library(argparser)
  library(ggplot2)
})
set.seed(1)
options(warn=1)

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

remove.covariate.effects.X <- function (X, Z) {
  # include the intercept term
  if (any(Z[,1]!=1)) Z = cbind(1, Z)
  A   <- forceSymmetric(crossprod(Z))
  SZX <- as.matrix(solve(A,t(Z) %*% X))
  X <- X - Z %*% SZX
  return(list(X=X, SZX=SZX))
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

p <- arg_parser("Run colocalization")
p <- add_argument(p, '--gwas_sumstatsFile', short='-s', help='Sumstats file with at least "beta" and "se" columns')
p <- add_argument(p, '--gwas_ldPrefix', short='-l', help='Prefix of GWAS LD matrix files.')
p <- add_argument(p, '--gwas_sampleSize', short='-n', help='GWAS sample size')
p <- add_argument(p, '--gwas_type', short='-t', help='GWAS type (quant or cc)')
p <- add_argument(p, '--eQTL_nominalResultsFile', short='-N', help='eQTL nominal results file from FastQTL')
p <- add_argument(p, '--eQTL_gtPrefix', short='-g', help='Prefix (including path) of MAGE genotype matrix files. There should be three files, ending ".012", ".012.indv", and ".012.snps"')
p <- add_argument(p, '--eQTL_covFile', short='-c', help='File with additional covariates to be regressed out')
p <- add_argument(p, '--outPrefix', short='-o', help="Prefix (including path) of output files")
p <- add_argument(p, '--overlap_vars', short='-v', help="Limit analysis to overlapping variants", flag=TRUE)
args <- parse_args(p)


#===============================#
# Read in and prepare GWAS data #
#===============================#

# Get GWAS sample size
gwas_sampleSize <- as.integer(args$gwas_sampleSize)

# Read in sumstats
cat('\t\t- Reading in GWAS sumstats\n')
gwas_sumstats <- read.table(args$gwas_sumstatsFile, sep='\t', header=TRUE)
rownames(gwas_sumstats) <- paste(gwas_sumstats$chr, gwas_sumstats$pos, gwas_sumstats$ref, gwas_sumstats$alt, sep='_')

# Get MAF
gwas_sumstats$raf <- (1-gwas_sumstats$af)
gwas_sumstats <- transform(gwas_sumstats, maf = pmin(raf, af))

# Read in LD matrix
cat('\t\t- Reading in LD matrix\n')
gwas_ldMatrixFile <- paste0(args$gwas_ldPrefix, '.ld')
gwas_ldVariantsFile <- paste0(args$gwas_ldPrefix, '.snps')
gwas_ld_matrix <- as.matrix(read.table(gwas_ldMatrixFile, sep='\t'))
gwas_ld_snps <- read.table(gwas_ldVariantsFile, sep='\t', header=TRUE)
rownames(gwas_ld_snps) <- paste(gwas_ld_snps$chr, gwas_ld_snps$pos, gwas_ld_snps$ref, gwas_ld_snps$alt, sep='_')
rownames(gwas_ld_matrix) <- rownames(gwas_ld_snps)
colnames(gwas_ld_matrix) <- rownames(gwas_ld_snps)

# Subset to shared variants (between sumstats and LD matrix)
gwas_snps <- rownames(gwas_ld_snps)[rownames(gwas_ld_snps) %in% rownames(gwas_sumstats)]
gwas_sumstats <- gwas_sumstats[gwas_snps, ]
gwas_ld_matrix <- gwas_ld_matrix[gwas_snps, gwas_snps]
gwas_ld_snps <- gwas_ld_snps[gwas_snps, ]

# Fill in missing values
if (sum(is.na(gwas_ld_matrix)) != 0) {
  cat('\t\t\t* Some LD values are missing. They have been set to 0\n')
  gwas_ld_matrix[is.na(gwas_ld_matrix)] <- 0
}


#===============================#
# Read in and prepare eQTL data #
#===============================#

# Read in nominal results
cat('\t\t- Reading in eQTL nominal results\n')
nominal_results <- read.table(args$eQTL_nominalResultsFile, sep='\t', header=TRUE)
nominal_results$variantFormattedID <- paste(nominal_results$variantChrom, nominal_results$variantPosition, nominal_results$variantRef, nominal_results$variantAlt, sep='_')

# Read in eQTL genotypes
cat('\t\t- Reading in eQTL genotypes\n')
eQTL_snps <- read.table(paste0(args$eQTL_gtPrefix, ".012.snps"))[,1]
sample_ids <- read.table(paste0(args$eQTL_gtPrefix, ".012.indv"))[,1]
X <- as.matrix(read.table(paste0(args$eQTL_gtPrefix, ".012"), sep='\t', header=FALSE, row.names=1))
rownames(X) <- sample_ids
colnames(X) <- eQTL_snps

# Get eQTL sample size
eQTL_sampleSize <- length(sample_ids)

# Read in covariates
cat('\t\t- Reading in covariates\n')
Z <- as.data.frame(t(read.table(args$eQTL_covFile, sep='\t', header=FALSE)))
colnames(Z) <- Z[1,]
rownames(Z) <- Z[,1]
Z <- Z[-1,-1]
Z <- Z[sample_ids,]
Z <- numeric.matrix(Z)

# Regress out covariates
cat('\t\t- Regressing out covariates\n')
out <- remove.covariate.effects.X(X, Z)
X <- out$X

# Calculate MAGE LD matrix
cat('\t\t- Calculating eQTL LD matrix\n')
eQTL_ld_matrix <- cor(X)


#==========================================#
# Subset to shared variants (if requested) #
#==========================================#

if (args$overlap_vars) {
  cat('\t\t- Subsetting to shared variants\n')
  eQTL_snps <- eQTL_snps[eQTL_snps %in% gwas_snps]
  nominal_results <- nominal_results[nominal_results$variantFormattedID %in% eQTL_snps, ]
  eQTL_ld_matrix <- eQTL_ld_matrix[eQTL_snps, eQTL_snps]

  gwas_snps <- gwas_snps[gwas_snps %in% eQTL_snps]
  gwas_sumstats <- gwas_sumstats[gwas_snps, ]
  gwas_ld_matrix <- gwas_ld_matrix[gwas_snps, gwas_snps]
}


#=================================#
# Check GWAS Z-score LD agreement #
#=================================#

# cat('\t\t- Checking agreement between Z-scores and LD for GWAS\n')

# gwas_sumstats$zscore <- gwas_sumstats$beta / gwas_sumstats$se

# # Get s-value
# s <- estimate_s_rss(gwas_sumstats$zscore, gwas_ld_matrix, n=gwas_sampleSize)
# cat(paste0('\t\t\t* lambda (lower better) = ', s, '\n'))

# # Get kriging plot
# cond_dist <- kriging_rss(gwas_sumstats$zscore, gwas_ld_matrix, n=gwas_sampleSize, s=s)
# plot <- cond_dist$plot

# outPlot <- paste0(args$outPrefix, ".kriging.png")
# ggsave(outPlot, plot=plot, width=5, height=5, units="in", dpi=500)


#===========================#
# Set up GWAS coloc dataset #
#===========================#

gwas_dataset <- list(beta = gwas_sumstats$beta,
                     varbeta = gwas_sumstats$se^2,
                     snp = gwas_snps,
                     position=gwas_sumstats$pos,
                     type=args$gwas_type,
                     N=gwas_sampleSize,
                     MAF=gwas_sumstats$maf,
                     LD=gwas_ld_matrix
                     )


#============================#
# Run SuSiE on GWAS sumstats #
#============================#

cat('\t\t- Running SuSiE on GWAS sumstats\n')
gwas_susie_results <- suppressMessages(suppressWarnings(runsusie(gwas_dataset, L=10, coverage=0.95, min_abs_corr=0.5, estimate_residual_variance=TRUE)))
numCS_GWAS <- length(gwas_susie_results$set$cs)


#=======================#
# Prepare output tables #
#=======================#

# Create full colocalization results output table
coloc_results_summary <- data.frame(ensemblID=character(),
                                    geneSymbol=character(),
                                    numCS_GWAS=numeric(),
                                    numCS_eQTL=numeric(),
                                    numColoc_0.5=numeric(),
                                    numColoc_0.8=numeric()
                                    )

coloc_results <- data.frame(eQTL_gene_ensemblID=character(),
                            eQTL_gene_geneSymbol=character(),
                            eQTL_variantChrom=character(),
                            eQTL_variantPosition=character(),
                            eQTL_variantRef=character(),
                            eQTL_variantAlt=character(),
                            eQTL_variant_kgpID=character(),
                            eQTL_variant_rsID=character(),
                            eQTL_geneCS=character(),
                            GWAS_variantChrom=character(),
                            GWAS_variantPosition=character(),
                            GWAS_variantRef=character(),
                            GWAS_variantAlt=character(),
                            GWAS_variant_kgpID=character(),
                            GWAS_variant_rsID=character(),
                            GWAS_CS=character(),
                            PP.H0=numeric(),
                            PP.H1=numeric(),
                            PP.H2=numeric(),
                            PP.H3=numeric(),
                            PP.H4=numeric()
                            )

eQTL_results <- data.frame(variantChrom=character(),
                           variantPosition=character(),
                           variantRef=character(),
                           variantAlt=character(),
                           variant_kgpID=character(),
                           variant_rsID=character(),
                           ensemblID=character(),
                           geneSymbol=character(),
                           variantPIP=numeric(),
                           variantCredibleSet=character()
                           )

GWAS_results <- data.frame(variantChrom=character(),
                           variantPosition=character(),
                           variantRef=character(),
                           variantAlt=character(),
                           variant_kgpID=character(),
                           variant_rsID=character(),
                           variantPIP=numeric(),
                           variantCredibleSet=character()
                           )


#===================================#
# Run colocalization for each eGene #
#===================================#

any_moderate_colocs <- FALSE

cat('\t\t- Running SuSiE on eQTL sumstats and doing colocalization\n')
for (ensemblID in unique(nominal_results$ensemblID)) {

  geneSymbol <- nominal_results[nominal_results$ensemblID == ensemblID, 'geneSymbol'][1]
  cat(paste0('\t\t\t* ', ensemblID, ' (', geneSymbol, ')\n'))
  gene_nominal_results <- nominal_results[which(nominal_results$ensemblID == ensemblID), ]
  rownames(gene_nominal_results) <- gene_nominal_results$variantFormattedID

  gene_ld_matrix <- eQTL_ld_matrix[gene_nominal_results$variantFormattedID, gene_nominal_results$variantFormattedID]

  gene_eqtl_dataset <- list(beta = gene_nominal_results$slope,
                            varbeta = gene_nominal_results$slope_se^2,
                            snp = gene_nominal_results$variantFormattedID,
                            position=gene_nominal_results$variantPosition,
                            type="quant",
                            N=eQTL_sampleSize,
                            MAF=gene_nominal_results$maf,
                            LD=gene_ld_matrix
                            )
  
  # Run SuSiE on eQTL results
  gene_eqtl_susie_results <- suppressMessages(suppressWarnings(runsusie(gene_eqtl_dataset, L=10, coverage=0.95, min_abs_corr=0.5, estimate_residual_variance=TRUE)))
  numCS_eQTL <- length(gene_eqtl_susie_results$set$cs)

  gene_coloc_res <- coloc.susie(gwas_susie_results, gene_eqtl_susie_results)
  gene_coloc_res_sum <- gene_coloc_res$summary

  # Add results to output tables
  if (is.null(gene_coloc_res_sum)) {
    numColoc_0.5 <- 0
    numColoc_0.8 <- 0
    coloc_results_summary[nrow(coloc_results_summary)+1, ] = c(ensemblID, geneSymbol, numCS_GWAS, numCS_eQTL, numColoc_0.5, numColoc_0.8)
    next
  }

  numColoc_0.5 <- sum(gene_coloc_res_sum$PP.H4.abf >= 0.5)
  numColoc_0.8 <- sum(gene_coloc_res_sum$PP.H4.abf >= 0.8)
  coloc_results_summary[nrow(coloc_results_summary)+1, ] <- c(ensemblID, geneSymbol, numCS_GWAS, numCS_eQTL, numColoc_0.5, numColoc_0.8)

  for (i in 1:nrow(gene_coloc_res_sum)) {
    GWAS_credibleSet <- paste0('L', gene_coloc_res_sum$idx1[i])
    GWAS_variantID <- gene_coloc_res_sum$hit1[i]
    GWAS_variant_kgpID <- gene_nominal_results[GWAS_variantID, 'variant_kgpID']
    GWAS_variant_rsID <- gene_nominal_results[GWAS_variantID, 'variant_rsID']
    GWAS_variantChrom <- gwas_sumstats[GWAS_variantID, 'chr']
    GWAS_variantPosition <- gwas_sumstats[GWAS_variantID, 'pos']
    GWAS_variantRef <- gwas_sumstats[GWAS_variantID, 'ref']
    GWAS_variantAlt <- gwas_sumstats[GWAS_variantID, 'alt']
    eQTL_credibleSet <- paste0('L', gene_coloc_res_sum$idx2[i])
    eQTL_variantID <- gene_coloc_res_sum$hit2[i]
    eQTL_variant_kgpID <- gene_nominal_results[eQTL_variantID, 'variant_kgpID']
    eQTL_variant_rsID <- gene_nominal_results[eQTL_variantID, 'variant_rsID']
    eQTL_variantChrom <- gene_nominal_results[eQTL_variantID, 'variantChrom']
    eQTL_variantPosition <- gene_nominal_results[eQTL_variantID, 'variantPosition']
    eQTL_variantRef <- gene_nominal_results[eQTL_variantID, 'variantRef']
    eQTL_variantAlt <- gene_nominal_results[eQTL_variantID, 'variantAlt']
    PP.H0 <- gene_coloc_res_sum$PP.H0.abf[i]
    PP.H1 <- gene_coloc_res_sum$PP.H1.abf[i]
    PP.H2 <- gene_coloc_res_sum$PP.H2.abf[i]
    PP.H3 <- gene_coloc_res_sum$PP.H3.abf[i]
    PP.H4 <- gene_coloc_res_sum$PP.H4.abf[i]

    coloc_results[nrow(coloc_results)+1, ] <- c(ensemblID, geneSymbol, eQTL_variantChrom, eQTL_variantPosition, eQTL_variantRef, eQTL_variantAlt, eQTL_variant_kgpID, eQTL_variant_rsID, eQTL_credibleSet, GWAS_variantChrom, GWAS_variantPosition, GWAS_variantRef, GWAS_variantAlt, GWAS_variant_kgpID, GWAS_variant_rsID, GWAS_credibleSet, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)
  }

  if (numColoc_0.5 > 0) {
    any_moderate_colocs <- TRUE

    for (i in 1:length(names(gene_eqtl_susie_results$set$cs))) {
      cred_set <- names(gene_eqtl_susie_results$set$cs)[i]
      for (var_index in gene_eqtl_susie_results$set$cs[cred_set][[1]]) {
        variantID <- gene_nominal_results$variantFormattedID[var_index]
        variant_kgpID <- gene_nominal_results[variantID, 'variant_kgpID']
        variant_rsID <- gene_nominal_results[variantID, 'variant_rsID']
        variantPIP <- gene_eqtl_susie_results$pip[var_index]
        variantChrom <- gene_nominal_results[variantID, 'variantChrom']
        variantPosition <- gene_nominal_results[variantID, 'variantPosition']
        variantRef <- gene_nominal_results[variantID, 'variantRef']
        variantAlt <- gene_nominal_results[variantID, 'variantAlt']
        eQTL_results[nrow(eQTL_results)+1, ] <- c(variantChrom, variantPosition, variantRef, variantAlt, variant_kgpID, variant_rsID, ensemblID, geneSymbol, variantPIP, cred_set)
      }
    }
  }
}


#==============#
# Write output #
#==============#

outColocSummaryFile <- paste0(args$outPrefix, ".coloc.summary.txt")
outColocResultsFile <- paste0(args$outPrefix, ".coloc.fullResults.txt")

write.table(coloc_results_summary, outColocSummaryFile, append=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, na="NA")
write.table(coloc_results, outColocResultsFile, append=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, na="NA")


#==================================================#
# Write credible sets for moderate colocalizations #
#==================================================#

if (any_moderate_colocs) {

  for (i in 1:length(names(gwas_susie_results$set$cs))) {
    cred_set <- names(gwas_susie_results$set$cs)[i]
    for (var_index in gwas_susie_results$set$cs[cred_set][[1]]) {
      variantID <- gwas_snps[var_index]
      variant_kgpID <- gene_nominal_results[variantID, 'variant_kgpID']
      variant_rsID <- gene_nominal_results[variantID, 'variant_rsID']
      variantPIP <- gwas_susie_results$pip[var_index]
      variantChrom <- gene_nominal_results[variantID, 'variantChrom']
      variantPosition <- gene_nominal_results[variantID, 'variantPosition']
      variantRef <- gene_nominal_results[variantID, 'variantRef']
      variantAlt <- gene_nominal_results[variantID, 'variantAlt']
      GWAS_results[nrow(GWAS_results)+1, ] <- c(variantChrom, variantPosition, variantRef, variantAlt, variant_kgpID, variant_rsID, variantPIP, cred_set)
    }
  }

  outGWAS_credSetFile <- paste0(args$outPrefix, ".coloc.GWAS_CS.txt")
  outeQTL_credSetFile <- paste0(args$outPrefix, ".coloc.eQTL_CS.txt")

  write.table(GWAS_results, outGWAS_credSetFile, append=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, na='NA')
  write.table(eQTL_results, outeQTL_credSetFile, append=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, na='NA')

}


#===========#
# Finish up #
#===========#

cat('\t\t- Done.\n')
