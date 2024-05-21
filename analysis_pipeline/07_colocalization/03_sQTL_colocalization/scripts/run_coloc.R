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
p <- add_argument(p, '--sQTL_nominalResultsFile', short='-N', help='sQTL nominal results file from FastQTL')
p <- add_argument(p, '--sQTL_gtPrefix', short='-g', help='Prefix (including path) of MAGE genotype matrix files. There should be three files, ending ".012", ".012.indv", and ".012.snps"')
p <- add_argument(p, '--sQTL_covFile', short='-c', help='File with additional covariates to be regressed out')
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
# Read in and prepare sQTL data #
#===============================#

# Read in nominal results
cat('\t\t- Reading in sQTL nominal results\n')
nominal_results <- read.table(args$sQTL_nominalResultsFile, sep='\t', header=TRUE)
nominal_results$variantFormattedID <- paste(nominal_results$variantChrom, nominal_results$variantPosition, nominal_results$variantRef, nominal_results$variantAlt, sep='_')

# Read in sQTL genotypes
cat('\t\t- Reading in sQTL genotypes\n')
sQTL_snps <- read.table(paste0(args$sQTL_gtPrefix, ".012.snps"))[,1]
sample_ids <- read.table(paste0(args$sQTL_gtPrefix, ".012.indv"))[,1]
X <- as.matrix(read.table(paste0(args$sQTL_gtPrefix, ".012"), sep='\t', header=FALSE, row.names=1))
rownames(X) <- sample_ids
colnames(X) <- sQTL_snps

# Get sQTL sample size
sQTL_sampleSize <- length(sample_ids)

# Read in covariates
cat('\t\t- Reading in covariates\n')
Z <- as.data.frame(t(read.table(args$sQTL_covFile, sep='\t', header=FALSE)))
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
cat('\t\t- Calculating sQTL LD matrix\n')
sQTL_ld_matrix <- cor(X)


#==========================================#
# Subset to shared variants (if requested) #
#==========================================#

if (args$overlap_vars) {
  cat('\t\t- Subsetting to shared variants\n')
  sQTL_snps <- sQTL_snps[sQTL_snps %in% gwas_snps]
  nominal_results <- nominal_results[nominal_results$variantFormattedID %in% sQTL_snps, ]
  sQTL_ld_matrix <- sQTL_ld_matrix[sQTL_snps, sQTL_snps]

  gwas_snps <- gwas_snps[gwas_snps %in% sQTL_snps]
  gwas_sumstats <- gwas_sumstats[gwas_snps, ]
  gwas_ld_matrix <- gwas_ld_matrix[gwas_snps, gwas_snps]
}


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
coloc_results_summary <- data.frame(intronID=character(),
                                    ensemblID=character(),
                                    geneSymbol=character(),
                                    numCS_GWAS=numeric(),
                                    numCS_sQTL=numeric(),
                                    numColoc_0.5=numeric(),
                                    numColoc_0.8=numeric()
                                    )

coloc_results <- data.frame(sQTL_intronID=character(),
                            sQTL_gene_ensemblID=character(),
                            sQTL_gene_geneSymbol=character(),
                            sQTL_variantChrom=character(),
                            sQTL_variantPosition=character(),
                            sQTL_variantRef=character(),
                            sQTL_variantAlt=character(),
                            sQTL_variant_kgpID=character(),
                            sQTL_variant_rsID=character(),
                            sQTL_geneCS=character(),
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

sQTL_results <- data.frame(variantChrom=character(),
                           variantPosition=character(),
                           variantRef=character(),
                           variantAlt=character(),
                           variant_kgpID=character(),
                           variant_rsID=character(),
                           intronID=character(),
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

cat('\t\t- Running SuSiE on sQTL sumstats and doing colocalization\n')
for (intronID in unique(nominal_results$intronID)) {

  ensemblID <- nominal_results[nominal_results$intronID == intronID, 'ensemblID'][1]
  geneSymbol <- nominal_results[nominal_results$intronID == intronID, 'geneSymbol'][1]
  cat(paste0('\t\t\t* ', intronID, ' (', geneSymbol, ')\n'))
  intron_nominal_results <- nominal_results[which(nominal_results$intronID == intronID), ]
  rownames(intron_nominal_results) <- intron_nominal_results$variantFormattedID

  intron_ld_matrix <- sQTL_ld_matrix[intron_nominal_results$variantFormattedID, intron_nominal_results$variantFormattedID]

  intron_eqtl_dataset <- list(beta = intron_nominal_results$slope,
                            varbeta = intron_nominal_results$slope_se^2,
                            snp = intron_nominal_results$variantFormattedID,
                            position=intron_nominal_results$variantPosition,
                            type="quant",
                            N=sQTL_sampleSize,
                            MAF=intron_nominal_results$maf,
                            LD=intron_ld_matrix
                            )
  
  # Run SuSiE on sQTL results
  intron_eqtl_susie_results <- suppressMessages(suppressWarnings(runsusie(intron_eqtl_dataset, L=10, coverage=0.95, min_abs_corr=0.5, estimate_residual_variance=TRUE)))
  numCS_sQTL <- length(intron_eqtl_susie_results$set$cs)

  intron_coloc_res <- coloc.susie(gwas_susie_results, intron_eqtl_susie_results)
  intron_coloc_res_sum <- intron_coloc_res$summary

  # Add results to output tables
  if (is.null(intron_coloc_res_sum)) {
    numColoc_0.5 <- 0
    numColoc_0.8 <- 0
    coloc_results_summary[nrow(coloc_results_summary)+1, ] = c(intronID, ensemblID, geneSymbol, numCS_GWAS, numCS_sQTL, numColoc_0.5, numColoc_0.8)
    next
  }

  numColoc_0.5 <- sum(intron_coloc_res_sum$PP.H4.abf >= 0.5)
  numColoc_0.8 <- sum(intron_coloc_res_sum$PP.H4.abf >= 0.8)
  coloc_results_summary[nrow(coloc_results_summary)+1, ] <- c(intronID, ensemblID, geneSymbol, numCS_GWAS, numCS_sQTL, numColoc_0.5, numColoc_0.8)

  for (i in 1:nrow(intron_coloc_res_sum)) {
    GWAS_credibleSet <- paste0('L', intron_coloc_res_sum$idx1[i])
    GWAS_variantID <- intron_coloc_res_sum$hit1[i]
    GWAS_variant_kgpID <- intron_nominal_results[GWAS_variantID, 'variant_kgpID']
    GWAS_variant_rsID <- intron_nominal_results[GWAS_variantID, 'variant_rsID']
    GWAS_variantChrom <- gwas_sumstats[GWAS_variantID, 'chr']
    GWAS_variantPosition <- gwas_sumstats[GWAS_variantID, 'pos']
    GWAS_variantRef <- gwas_sumstats[GWAS_variantID, 'ref']
    GWAS_variantAlt <- gwas_sumstats[GWAS_variantID, 'alt']
    sQTL_credibleSet <- paste0('L', intron_coloc_res_sum$idx2[i])
    sQTL_variantID <- intron_coloc_res_sum$hit2[i]
    sQTL_variant_kgpID <- intron_nominal_results[sQTL_variantID, 'variant_kgpID']
    sQTL_variant_rsID <- intron_nominal_results[sQTL_variantID, 'variant_rsID']
    sQTL_variantChrom <- intron_nominal_results[sQTL_variantID, 'variantChrom']
    sQTL_variantPosition <- intron_nominal_results[sQTL_variantID, 'variantPosition']
    sQTL_variantRef <- intron_nominal_results[sQTL_variantID, 'variantRef']
    sQTL_variantAlt <- intron_nominal_results[sQTL_variantID, 'variantAlt']
    PP.H0 <- intron_coloc_res_sum$PP.H0.abf[i]
    PP.H1 <- intron_coloc_res_sum$PP.H1.abf[i]
    PP.H2 <- intron_coloc_res_sum$PP.H2.abf[i]
    PP.H3 <- intron_coloc_res_sum$PP.H3.abf[i]
    PP.H4 <- intron_coloc_res_sum$PP.H4.abf[i]

    coloc_results[nrow(coloc_results)+1, ] <- c(intronID, ensemblID, geneSymbol, sQTL_variantChrom, sQTL_variantPosition, sQTL_variantRef, sQTL_variantAlt, sQTL_variant_kgpID, sQTL_variant_rsID, sQTL_credibleSet, GWAS_variantChrom, GWAS_variantPosition, GWAS_variantRef, GWAS_variantAlt, GWAS_variant_kgpID, GWAS_variant_rsID, GWAS_credibleSet, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)
  }

  if (numColoc_0.5 > 0) {
    any_moderate_colocs <- TRUE

    for (i in 1:length(names(intron_eqtl_susie_results$set$cs))) {
      cred_set <- names(intron_eqtl_susie_results$set$cs)[i]
      for (var_index in intron_eqtl_susie_results$set$cs[cred_set][[1]]) {
        variantID <- intron_nominal_results$variantFormattedID[var_index]
        variant_kgpID <- intron_nominal_results[variantID, 'variant_kgpID']
        variant_rsID <- intron_nominal_results[variantID, 'variant_rsID']
        variantPIP <- intron_eqtl_susie_results$pip[var_index]
        variantChrom <- intron_nominal_results[variantID, 'variantChrom']
        variantPosition <- intron_nominal_results[variantID, 'variantPosition']
        variantRef <- intron_nominal_results[variantID, 'variantRef']
        variantAlt <- intron_nominal_results[variantID, 'variantAlt']
        sQTL_results[nrow(sQTL_results)+1, ] <- c(variantChrom, variantPosition, variantRef, variantAlt, variant_kgpID, variant_rsID, intronID, ensemblID, geneSymbol, variantPIP, cred_set)
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
      variant_kgpID <- intron_nominal_results[variantID, 'variant_kgpID']
      variant_rsID <- intron_nominal_results[variantID, 'variant_rsID']
      variantPIP <- gwas_susie_results$pip[var_index]
      variantChrom <- intron_nominal_results[variantID, 'variantChrom']
      variantPosition <- intron_nominal_results[variantID, 'variantPosition']
      variantRef <- intron_nominal_results[variantID, 'variantRef']
      variantAlt <- intron_nominal_results[variantID, 'variantAlt']
      GWAS_results[nrow(GWAS_results)+1, ] <- c(variantChrom, variantPosition, variantRef, variantAlt, variant_kgpID, variant_rsID, variantPIP, cred_set)
    }
  }

  outGWAS_credSetFile <- paste0(args$outPrefix, ".coloc.GWAS_CS.txt")
  outsQTL_credSetFile <- paste0(args$outPrefix, ".coloc.sQTL_CS.txt")

  write.table(GWAS_results, outGWAS_credSetFile, append=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, na='NA')
  write.table(sQTL_results, outsQTL_credSetFile, append=FALSE, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, na='NA')

}


#===========#
# Finish up #
#===========#

cat('\t\t- Done.\n')
