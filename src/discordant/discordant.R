library(Biobase)
library(discordant)
library(tidyverse)
library(data.table)
library(optparse)


# parsing arguments
opts <- list(
    make_option(c("-a", "--input1"),
        type = "character",
        help = "Input file1"
    ),
    make_option(c("-b", "--input2"),
        type = "character",
        help = "Input file2"
    ),
    make_option(c("-o", "--output_path"),
        type = "character",
        help = "Output path"
    ),
    make_option(c("-c", "--correlation_method"),
        type = "character",
        default = "spearman",
        help = "Correlation method to be used (spearman, pearson, bwmc, sparcc)"
    ),
    make_option(c("-n", "--components"),
        type = "numeric",
        default = 3,
        help = "number of mixture components to be fitted to correlation coefficients"
    ),
    make_option(c("-s", "--subsampling"),
        type = "logical",
        default = TRUE,
        help = "use subsampling during fitting with EM algorithm, recommended for large datasets."
    ),
    make_option(c("-t", "--transform"),
        type = "logical",
        default = TRUE,
        help = "Fisher z-transform correlation coefficients before fitting mixture model."
    ),
    make_option(c("-p", "--pp_threshold"),
        type = "numeric",
        default = 0.95,
        help = "Fisher z-transform correlation coefficients before fitting mixture model."
    ),
    make_option(c("-w", "--abs_weight_thresh"),
        type = "numeric",
        default = 0.6,
        help = "Filter by absolute weight (correlation difference across the conditions) threshold."
    )
)

opts <- parse_args(OptionParser(option_list = opts))


# getting condition names --------
cond_name_1 <- sub(".*out_", "", opts$input1)
cond_name_1 <- sub(".tsv", "", cond_name_1)

cond_name_2 <- sub(".*out_", "", opts$input2)
cond_name_2 <- sub(".tsv", "", cond_name_2)


# reading data
data_1 <- fread(opts$input1, data.table = FALSE)
data_2 <- fread(opts$input2, data.table = FALSE)


# creating groups based on data --
groups <- c(rep(1, times = ncol(data_1) - 1), rep(2, times = ncol(data_2) - 1)) # -1 to exclude Gene column


# merging data
data <- bind_cols(data_1, select(data_2, -Gene))
rownames(data) <- data[, 1]


# creating ExpressionSet
data <- data %>%
    select(-Gene) %>%
    as.matrix() %>%
    ExpressionSet() # ExpressionSet format required for discordant


# calculating correlation vectors
correlation_vectors <- createVectors(
    x = data,
    groups = groups,
    cor.method = opts$correlation_method # "spearman", "pearson", "bwmc", "sparcc"
)


# calculating correlation difference to be considered weight
corrdiff <- correlation_vectors$v1 - correlation_vectors$v2


# running discordant algorithm
results <- discordantRun(
    v1 = correlation_vectors$v1,
    v2 = correlation_vectors$v2,
    x = data,
    components = opts$components, # number of mixture components to be fitted to correlation coefficients
    subsampling = opts$subsampling, # use subsampling during fitting with EM algorithm, recommended for large datasets.
    iter = 100, # number of iterations for subsampling
    transform = opts$transform # Fisher z-transform correlation vectors
)


# processing and saving results
if (length(results$discordPPVector) > 0) {
    df <- data.frame(names = names(results$discordPPVector), pp = results$discordPPVector, weight = corrdiff)
    df <- df %>%
        separate(names, into = c("target", "regulator"), sep = "_") %>%
        filter(pp >= opts$pp_threshold) %>% # filtering by posterior probability
        filter(abs(weight) >= opts$abs_weight_thresh) %>% # filtering by correlation difference
        mutate(condition = if_else(weight >= 0, cond_name_1, cond_name_2)) %>% # assigning condition based on weight(correlation difference)
        select(target, regulator, condition, weight)
    fwrite(df, paste0(opts$output_path, "/network.tsv"), sep = "\t")
} else {
    print("No discordant genes found")
    df <- data.frame(target = NULL, source = NULL, weight = NULL, condition = NULL)
    fwrite(df, paste0(opts$output_path, "/network.tsv"), sep = "\t")
}
