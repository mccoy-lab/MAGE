#!/usr/bin/env Rscript

suppressMessages(library(leafcutter))
library(argparser)

#===================#
# Collect arguments #
#===================#

p <- arg_parser("LeafCutter: map clusters to genes")
p <- add_argument(p, "leafcutterQuantsFile", help="Intron counts/ratios file from LeafCutter, typically <prefix>_perind_numers.counts.gz or <prefix>_perind.counts.gz, respectively")
p <- add_argument(p, "exonFile", help="File listing all unique exons in annotation. Must have columns: chr, start, end, strand, gene_id[, gene_name].")
p <- add_argument(p, "outFile", help="Output file name")
args <- parse_args(p)


#==============#
# Read in data #
#==============#

cat("LeafCutter: mapping clusters to genes\n")
intronQuants <- read.table(args$leafcutterQuantsFile, header=TRUE, check.names=FALSE, row.names=1)
intronMeta <- get_intron_meta(rownames(intronQuants))

exonTable <- read.table(args$exonFile, header=TRUE, stringsAsFactors=FALSE)
stopifnot(is.element('gene_id', colnames(exonTable)))
exonTable[, 'gene_name'] <- exonTable[, 'gene_id']


#=======================#
# Map clusters to genes #
#=======================#

m <- map_clusters_to_genes(intronMeta, exonTable)
write.table(m, args$outFile, sep = "\t", quote=FALSE, row.names=FALSE)
