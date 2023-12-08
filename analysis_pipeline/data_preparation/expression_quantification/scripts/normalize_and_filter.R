#!/usr/bin/env Rscript

library(tximport)
library(edgeR)
library(matrixStats)
library(argparser)
suppressMessages(library(DESeq2))
library(BiocParallel)

#===================#
# Collect arguments #
#===================#

p <- arg_parser("Collate salmon expression quantifications, filter based on expression thresholds, and tmm normalize")
p <- add_argument(p, '--salmonFiles', short='-s', help='A file with two columns: 1) sampleID and 2) the filepath of the corresponding salmon "quant.sf" file')
p <- add_argument(p, '--tx2gene', short='-t', help='Comma separated file containing transcript to gene mappings for tximport')
p <- add_argument(p, '--keepGenes', short='-k', help='Column file with list of genes to limit output to (e.g. so as to filter out genes on chrY)')
p <- add_argument(p, '--outPrefix', short='-o', help="Prefix (with output directory) of output files")
args <- parse_args(p)


#====================#
# Read in expression #
#====================#

filesDF <- read.table(args$salmonFiles, header=FALSE, sep='\t')
files <- filesDF[,2]
names(files) <- filesDF[,1]

tx2gene <- read.csv(args$tx2gene)

txi <- tximport(files, type='salmon', tx2gene=tx2gene, abundanceCol='TPM', countsCol='NumReads', lengthCol='Length')


#==================================#
# Filter out lowly expressed genes #
#==================================#

# filter to genes with counts >= 6 and TPM >= 0.1 in at least 20% of samples
cts <- txi$counts
abd <- txi$abundance

keep_cts <- rowCounts(cts >= 6, value=TRUE) >= 0.2*ncol(cts)
keep_abd <- rowCounts(abd >= 0.1, value=TRUE) >= 0.2*ncol(cts)
keep_both <- keep_cts & keep_abd


#=============================#
# Convert pseudocounts to TMM #
#=============================#

# TMM normalize
y <- DGEList(cts)
y <- calcNormFactors(y, method='TMM')
tmm <- cpm(y)

# DESeq2 Normalize
dds <- DESeqDataSetFromTximport(txi, colData=as.data.frame(files), design=~1)
register(MulticoreParam(4))
dds <- suppressWarnings(DESeq(dds, parallel=TRUE))
normed_cts <- counts(dds, normalized=TRUE)
logged_normed_cts <- log2(normed_cts + 1)

# Apply mask
filtered_cts <- cts[keep_both, ]
filtered_tmm <- tmm[keep_both, ]
filtered_logged_normed_cts <- logged_normed_cts[keep_both, ]

# Limit to chosen genes
if (!is.null(args$keepGenes)) {
	keepGenes <- scan(args$keepGenes, what="character")
	filtered_cts <- filtered_cts[rownames(filtered_cts) %in% keepGenes, ]
	filtered_tmm <- filtered_tmm[rownames(filtered_tmm) %in% keepGenes, ]
	filtered_logged_normed_cts <- filtered_logged_normed_cts[rownames(filtered_logged_normed_cts) %in% keepGenes, ]
}

# Output
outPseudocountsFile <- paste(args$outPrefix, ".filtered.pseudocounts.tab", sep='')
outTMMFile <- paste(args$outPrefix, ".filtered.tmm.tab", sep='')
outCorrectedCountsFile <- paste(args$outPrefix, ".filtered.DESeq2-normalized-logged-pseudocounts.tab", sep='')

write.table(filtered_cts, file=outPseudocountsFile, quote=FALSE, sep='\t', col.names=NA)
write.table(filtered_tmm, file=outTMMFile, quote=FALSE, sep='\t', col.names=NA)
write.table(filtered_logged_normed_cts, file=outCorrectedCountsFile, quote=FALSE, sep='\t', col.names=NA)
