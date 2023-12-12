#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(manta)
	library(argparser)
	library(doParallel)
	library(itertools)
})


#==================#
# Define functions #
#==================#

# The default manta function hits errors when calculating small p-values.
# This version of the function is identical to the official manta function
# but does not calculate p-values
manta <- function(formula, data, transform = "none", type = "II", 
								contrasts = NULL, subset = NULL, fit = FALSE){
	
	## Checks
	transform <- match.arg(transform, c("none", "sqrt", "log"))
	type <- match.arg(type, c("I", "II", "III"))
	
	## Save call, build model frame, obtain responses
	cl <- match.call()
	m <- match(c("formula", "data"), names(cl), 0L)
	mf <- cl[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf$na.action = "na.omit"
		# The rows with at least one NA either in Y or X 
		# (only considering variables used in the formula) 
		# will be removed before transforming/centering
	mf[[1L]] <- quote(stats::model.frame)
	mf <- eval(mf, parent.frame())               
	mt <- attr(mf, "terms")
	response <- model.response(mf, "numeric")
	
	## Checks
	if (NCOL(response) < 2){
		stop("The number of response variables should be >= 2")
	}
	if (length(attr(mt, "term.labels")) < 1) {
		stop("The model should contain at least one predictor (excluding the intercept)")
	}
	
	## Transform and center responses, update model frame
	if(transform == "none"){
		Y <- response
	} else if (transform == "sqrt"){
		if (any(response < 0)) {
			stop("'sqrt' transformation requires all response values >= 0")
		}
		Y <- sqrt(response)
	} else if (transform == "log"){
		if (any(response <= 0)) {
			stop("'log' transformation requires all response values > 0")
		}
		Y <- log(response)
	}
	Y <- scale(Y, center = TRUE, scale = FALSE)
	mf[[1L]] <- Y
	
	## Define contrasts
	if(is.null(contrasts)){
		contrasts <- list(unordered = "contr.sum", ordered = "contr.poly")
		dc <- attr(mt, "dataClasses")[-1]
		contr.list <- lapply(dc, FUN = function(k){
			# No contrast for quantitative predictors
			# Sum contrasts for unordered categorical predictors
			# Polynomial contrasts for ordered categorical predictors
			contr.type <- switch(k, "factor" = contrasts$unordered,
													 "ordered" = contrasts$ordered)
			return(contr.type)
		})
		contr.list <- contr.list[!unlist(lapply(contr.list, is.null))]
	} else {
		contr.list <- contrasts
	}

	## Build model matrix
	X <- model.matrix(mt, mf, contr.list)
	
	## Fit lm
	lmfit <- lm.fit(X, Y)
	class(lmfit) <- c("mlm", "lm")
	lmfit$na.action <- attr(mf, "na.action")
	lmfit$contrasts <- attr(X, "contrasts")
	lmfit$xlevels <- .getXlevels(mt, mf)
	lmfit$call <- cl[c(1L, m)]
	lmfit$call[[1L]] <- quote(lm)
	if(length(contr.list) > 0) lmfit$call$contrasts <- quote(contr.list)
	lmfit$terms <- mt
	lmfit$model <- mf

	## Compute sums of squares, df's, pseudo-F statistics, partial R2s and eigenvalues 
	stats <- manta:::manta.ss(fit = lmfit, X = X, type = type, subset = subset)
	SS <- stats$SS
	df <- stats$df
	f.tilde <- stats$f.tilde
	r2 <- stats$r2
	e <- stats$e

	## Compute P-values
	l <- length(df) # SS[l], df[l] correspond to Residuals 
	# pv.acc <- mapply(manta:::p.asympt, ss = SS[-l], df = df[-l], MoreArgs = list(lambda = e))  
	pv.acc <- tryCatch({
		mapply(manta:::p.asympt, ss = SS[-l], df = df[-l], MoreArgs = list(lambda = e))
		}, error = function(e) {
			NAmatrix <- matrix(data=NA, nrow=2, ncol=l-1)
			colnames(NAmatrix) <- names(df)[1:l-1]
			return(NAmatrix)
		})

	## ANOVA table
	stats.l <- list(df, SS, SS/df, f.tilde, r2, pv.acc[1, ])
	# stats.l <- list(df, SS, SS/df, f.tilde, r2)
	cmat <- data.frame()
	for(i in seq(along = stats.l)) {
		for(j in names(stats.l[[i]])){
			cmat[j, i] <- stats.l[[i]][j]
		}
	}
	cmat <- as.matrix(cmat)
	colnames(cmat) <- c("Df", "Sum Sq", "Mean Sq", "F value", "R2", "Pr(>F)")

	## Output
	out <- list("call" = cl,
							"aov.tab" = cmat,
							"type" = type,
							# "precision" = pv.acc[2, ],
							"transform" = transform,
							"na.omit" = lmfit$na.action)
	if(fit){
		out$fit <- lmfit
	}

	## Update class
	class(out) <- c('manta', class(out))
	return(out)
}

run_manta <- function(usage_df) {
	
	PVE_df <- data.frame(matrix(ncol=7, nrow=0, dimnames=list(NULL, c("clusterID", "total_SS", "batch_SS", "sex_SS", "batch_sex_residual_SS", "continentalGroup_SS", "population_SS"))))

	for (target_cluster in unique(usage_df$cluster)) {

		# Run regression on each cluster
		target_cluster_df <- usage_df[which(usage_df['cluster'] == target_cluster), ]
		target_cluster_df_t <- t(target_cluster_df[, 3:ncol(target_cluster_df)])
		colnames(target_cluster_df_t) <- as.vector(target_cluster_df[, 'intron.cluster'])

		# Convert fractions to proportions
		proportion_cluster_df <- apply(target_cluster_df_t, c(1,2), function(x) eval(parse(text = x)))
		
		# Remove samples with 0 counts
		prepared_cluster_df <- na.omit(proportion_cluster_df)
		cluster_cov_df <- tx_cov_df[rownames(prepared_cluster_df), ]

		# Step 1: Regress out batch and sex
		confounding_results <- manta(prepared_cluster_df ~ ., data = subset(cluster_cov_df, select=c('batch', 'sex')), fit=TRUE, transform='sqrt')
		confounding_aov <- confounding_results$aov.tab

		batch_SS <- confounding_aov['batch', 'Sum Sq']
		sex_SS <- confounding_aov['sex', 'Sum Sq']
		residual_SS <- confounding_aov['Residuals', 'Sum Sq']
		total_SS <- batch_SS + sex_SS + residual_SS

		confounding_residuals <- confounding_results$fit$residual

		# Step 2: Run regression on residuals from Step 1 with continentalGroup
		continentalGroup_results <- manta(confounding_residuals ~ ., data = subset(cluster_cov_df, select='continentalGroup'), fit=TRUE)
		continentalGroup_aov <- continentalGroup_results$aov.tab

		continentalGroup_SS <- continentalGroup_aov['continentalGroup', 'Sum Sq']

		# Step 3: Run regression on residuals from Step 1 with population
		population_results <- manta(confounding_residuals ~ ., data=subset(cluster_cov_df, select='population'), fit=TRUE)
		population_aov <- population_results$aov.tab

		population_SS <- population_aov['population', 'Sum Sq']

		# Add cluster to PVE dataframe
		PVE_df[target_cluster, ] = c(target_cluster, total_SS, batch_SS, sex_SS, residual_SS, continentalGroup_SS, population_SS)

	}

	return(PVE_df)
}


#===================#
# Collect arguments #
#===================#

p <- arg_parser("Calculate proportion of splicing varition explained by population")
p <- add_argument(p, "--splicing_frac", short="-f", help="Intron exicision fractions from Leafcutter (the unfiltered file is typically \"[prefix]_perind.counts.gz\"). These are expected to be filtered prior to running this script. No further filtering is done in this script.")
p <- add_argument(p, "--covariates", short="-c", help="Table of sample covariates (batch, sex, continentalGroup, population). This should be for the same set of sample in the leafcutter counts file.")
p <- add_argument(p, "--threads", short="-t", type="integer", default=1, help="Number of threads to run on")
p <- add_argument(p, "--output", short="-o", help="Name of file to write output table to. This table will have a breakdown of sum of squares for each cluster.")
args <- parse_args(p)


#===============#
# Load in usage #
#===============#

cat('\n\tReading in data... ')

tx_usage_df <- read.table(file=args$splicing_frac, sep='\t', header=TRUE)

clusters <- unique(tx_usage_df$cluster)


#====================#
# Load in covariates #
#====================#

tx_cov_df <- read.table(file=args$covariates, sep='\t', header=TRUE, row.names=1)

tx_cov_df$batch <- as.factor(tx_cov_df$batch)
tx_cov_df$sex <- as.factor(tx_cov_df$sex)
tx_cov_df$continentalGroup <- as.factor(tx_cov_df$continentalGroup)
tx_cov_df$population <- as.factor(tx_cov_df$population)

tx_cov_df <- tx_cov_df[colnames(tx_usage_df)[3:ncol(tx_usage_df)], ]

cat('Done.\n')


#===========#
# Run ANOVA #
#===========#

cat('\n\tRunning MANTA for ', length(clusters), ' clusters (', nrow(tx_usage_df), ' introns) on ', args$threads, ' threads... ', sep='')

start.time <- Sys.time()

# Create chunks
chunks <- foreach (clusters=isplitVector(clusters, chunks=args$threads)) %do% {
	tx_usage_df[tx_usage_df$cluster %in% clusters, ]

}

# Remove full data.frame from memory
rm(tx_usage_df)

# Run chunks in parallel
cl <- makeCluster(args$threads)
registerDoParallel(cl)

PVE_df <- foreach(i=1:args$threads, .combine=rbind) %dopar% {
	run_manta(chunks[[i]])
}

stopCluster(cl)

end.time <- Sys.time()
tot.time <- round(difftime(end.time, start.time, units='secs'), 2)

cat('Done.\n\tTotal time: ', tot.time, 's\n', sep='')
cat('\n\tWriting output... ')

# Write to output
write.table(PVE_df, file=args$output, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

cat('Done.\n\n\tAll done!\n\n')
