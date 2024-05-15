library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(optparse)

# Load Network
opts <- list(
    make_option(c("-n", "--network"),
        type = "character",
        help = "Network file"
    ),
    make_option(c("-o", "--output_path"),
        type = "character",
        help = "Output path"
    )
)

opts <- parse_args(OptionParser(option_list = opts))

df <- fread(opts$network, data.table = FALSE)
df_conds <- split(df, df$condition)

genes <- unique(c(df$target, df$regulator))
genes_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)


# GO enrichment analysis

res_go <- lapply(df_conds, function(x) {
    unique_genes <- unique(c(x$target, x$regulator))

    ego <- enrichGO(
        gene = genes_entrez[genes_entrez$SYMBOL %in% unique_genes, ]$ENTREZID,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP", # Biological Processes
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
    )
})
res_go <- lapply(res_go, as.data.frame)
res_go <- do.call(rbind, res_go)
res_go <- res_go %>%
    mutate(condition = sub("\\.GO.*", "", rownames(res_go)))

if (nrow(res_go) > 0) {
    fwrite(res_go, paste0(opts$output_path, "/GO_enriched_terms.tsv"), sep = "\t")
}


# KEGG enrichment analysis

res_kegg <- lapply(df_conds, function(x) {
    unique_genes <- unique(c(x$target, x$regulator))

    kegg_results <- enrichKEGG(
        gene = genes_entrez[genes_entrez$SYMBOL %in% unique_genes, ]$ENTREZID,
        organism = "hsa",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
    )
})
res_kegg <- lapply(res_kegg, as.data.frame)
res_kegg <- do.call(rbind, res_kegg)
res_kegg <- res_kegg %>%
    mutate(condition = sub("\\.hsa.*", "", rownames(res_kegg)))

if (nrow(res_kegg) > 0) {
    fwrite(res_kegg, paste0(opts$output_path, "/KEGG_enriched_pathways.tsv"), sep = "\t")
}
