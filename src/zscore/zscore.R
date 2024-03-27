#!/usr/bin/env Rscript

library(optparse)
library(data.table)
library(igraph)
library(dcanr)

# parsing arguments --------------
option_list <- list(
    make_option(c("-c", "--input.file.1"), type = 'character', 
                help='Input file'),
    make_option(c("-d", "--input.file.2"), type = 'character',
                help='Input file'),
    make_option(c("-o", "--output.path"), type = 'character',
                help='output path')
)

opt <- parse_args(OptionParser(option_list=option_list))

# getting condition names --------
cond_name_1 <- sub('.*out_', "", opt$input.file.1)
cond_name_2 <- sub('.*out_', "", opt$input.file.2)

cond_name_1 <- sub(".tsv", "", cond_name_1)
cond_name_2 <- sub(".tsv", "", cond_name_2)

data.1 <- read.table(file = opt$input.file.1, sep = "\t", header = TRUE)
data.2 <- read.table(file = opt$input.file.2, sep = "\t", header = TRUE)

conditions <- c(rep(1, times=ncol(data.1)-1), rep(2, times=ncol(data.2)-1)) # -1 because of GENE column that has to be ignored

# merging data -------------------
data <- cbind(data.1, data.2[, 2:ncol(data.2)])
rownames(data) <- data[, 1]
data <- data[, 2:ncol(data)]
names(conditions) <- colnames(data)

# running diffcoex ---------------
scores <- dcScore(data, conditions, dc.method = 'zscore', cor.method = 'spearman')
raw_p <- dcTest(scores, data, conditions)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')
dcnet <- dcNetwork(scores, adj_p)
edgedf <- as_data_frame(dcnet, what = 'edges')

# catching edge case of empty df - 
if (nrow(edgedf) > 0) {
    edgedf$condition <- ifelse(edgedf$score >= 0, cond_name_1, cond_name_2)   

    edgedf <- edgedf[, c(1,2,5,3)]
    colnames(edgedf) <- c("target", "regulator", "condition", "weight")

    fwrite(edgedf, paste0(opt$output.path, "/network.txt"), sep="\t")
} else {
    edgedf <- data.frame(target = NULL, source = NULL, weight = NULL, condition = NULL)
    fwrite(edgedf, paste0(opt$output.path, "/network.txt"), sep="\t")
}