library(data.table)
library(tidyverse)
library(optparse)

source('12918_2013_1140_MOESM1_ESM.r')


# parsing arguments --------------
option_list <- list(
  make_option(c("-c", "--input.file.1"), type = 'character',
              help = 'Input file'),
  make_option(c("-d", "--input.file.2"), type = 'character',
              help = 'Input file'),
  make_option(c("-o", "--output.path"), type = 'character',
              help = 'output path'),
  make_option(
    "--dthresh",
    type = 'numeric',
    default = 0.3,
    help = 'Local topology analysis: threshold for topological dissimilarity'
  ),
  make_option(
    "--pthresh",
    type = 'numeric',
    default = 0.05,
    help = 'Global topology analysis: threshold for p-value'
  ),
  make_option(
    "--n.samples",
    type = 'numeric',
    default = 100,
    help = 'Number of samples for repeated subsampling'
  ),
  make_option(
    "--n.repetitions",
    type = 'numeric',
    default = 100,
    help = 'Number of repetitions for repeated subsampling'
  )
)


opt <- parse_args(OptionParser(option_list = option_list))






#### Functions -------------------------------------

# Randomly subsample to generate the null distribution

generate.random.subsamples <- function(data, n.samples) {
  # data: data to subsample
  # n.samples: subsampled size
  subsample.indices <-
    sample(ncol(data), size = n.samples, replace = FALSE)
  subsample.data <- data[, ..subsample.indices]
  return(subsample.data)
}


# repeat subsampling and detection of diff correlated genes for n.repetitions
perform.repetitions.local <-
  function(combined.data,
           n.repetitions,
           n.samples,
           path,
           file.name) {
    # combined.data: data of both conditions
    # n.repetitions: number of times subsampling and evaluation is repeated
    # n.samples: number of samples
    # path: output path to store DCloc files
    # file.name: file name of to save DCloc file
    
    for (i in 1:n.repetitions) {
      # Generate null distribution for a single repetition
      Mat.A <-
        generate.random.subsamples(combined.data[, -1], n.samples)
      Mat.B <-
        generate.random.subsamples(combined.data[, -1], n.samples)
      
     
      # Compare differential genes for a single repetition
      dcloc <- DCloc(Mat.A, Mat.B)
      write.table(
        dcloc,
        file = paste(path, file.name, i, ".tsv", sep = ""),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      
    }
    
  }

perform.repetitions.global <-
  function(combined.data,
           n.repetitions,
           n.samples,
           path,
           file.name) {
    # combined.data: data of both conditions
    # n.repetitions: number of times subsampling and evaluation is repeated
    # n.samples: number of samples
    # path: output path to store DCglob files
    # file.name: file name of to save DCglob file
    
    for (i in 1:n.repetitions) {
      Mat.A <- generate.random.subsamples(combined.data[, -1], n.samples)
      Mat.B <-
        generate.random.subsamples(combined.data[, -1], n.samples)
      
      # compute differential genes for subsets
      dcglob <- DCglob(Mat.A, Mat.B)
      write.table(
        dcglob,
        file = paste(path, file.name, i, ".tsv", sep = ""),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      
    }
    
  }


# calculate number of differentially correlated genes of local topology analysis
local.nr.genes <- function(df, d.threshold = 0.25) {
  # Parameters:
  # df: dataframe of local topology analysis; DCloc output
  # d.threshold: threshold to choose cutoff of differentially correlated genes based on topological dissimilarity
  
  # filter genes important for given condition
  df$keepGene.d <- (df$abs.top.dissim > d.threshold)
  
  return(sum(df$keepGene.d))
  
}


# aggregeate the results of the individual repetitions to estimate n0
local.n0.estimation <-
  function(n.repetitions,
           path,
           file.name,
           d.threshold = 0.25) {
    # n.repetitions: number of repetitions subsampling was conducted
    # path: output path where DCloc files are saved
    # d.threshold: cutoff threshold
    
    n0.results <- numeric(n.repetitions)
    
    
    for (i in 1:n.repetitions) {
      df <-
        fread(file = paste(path, file.name, i, ".tsv", sep = ""),
              sep = "\t")
      n0.results[i] <- local.nr.genes(df, d.threshold)
      
    }
    
    # return n0 estimation for both conditions
    return(n0.results)
    
  }



global.nr.genes <- function(df, p.threshold = 0.1) {
  # Parameters:
  # df: dataframe of global topology analysis, DCglob output
  # condition: condition of which the differentially correlated genes should be analyzed
  # p.threshold: threshold to choose cutoff of differentially correlated genes
  
  # filter genes important for given condition
  
  df$keepGene.p <- (df$`p-value` < p.threshold)
  return(sum(df$keepGene.p))
  
}


global.n0.estimation <-
  function(n.repetitions,
           path,
           file.name,
           p.threshold = 0.1) {
    # n.repetitions: number of repetitions subsampling was conducted
    # path: output path where DCloc files are saved
    # p.threshold: cutoff threshold
    
    n0.results <- numeric(n.repetitions)
    
    for (i in 1:n.repetitions) {
      df <-
        fread(file = paste(path, file.name, i, ".tsv", sep = ""),
              sep = "\t")
      n0.results[i] <- global.nr.genes(df, p.threshold)
    }
    
    return(n0.results)
    
  }


estimate.null.expectation <- function(n0.results) {
  return(list(
    mean = mean(n0.results),
    se = sd(n0.results) / sqrt(length(n0.results))
  ))
}


calculate.FDR <- function(nAB, n0, n.genes) {
  # FDR = nAB / ~n0
  return(nAB / (n.genes - n0))
}


cond.1 <- fread(file = opt$input.file.1,
                header = TRUE,
                sep = "\t")
cond.2 <- fread(file = opt$input.file.2,
                header = TRUE,
                sep = "\t")


# combine dataset for FDR

combined <- cbind(cond.1[, -1], cond.2[, -1])
shuffle.idx <-  sample(ncol(combined))
combined <- combined[, ..shuffle.idx]
combined <- cbind(cond.1$Gene, combined)
colnames(combined)[1] <- 'Gene'

n.repetitions <- opt$n.repetitions
n.samples <- opt$n.samples

data.cond.1 <- generate.random.subsamples(cond.1[, -1], n.samples)
data.cond.2 <- generate.random.subsamples(cond.2[, -1], n.samples)

dcloc.AB <- DCloc(data.cond.1, data.cond.2)
dcglob.AB <- DCglob(data.cond.1, data.cond.2)


n.repetitions <- opt$n.repetitions
n.samples <- opt$n.samples

data.cond.1 <- generate.random.subsamples(cond.1[, -1], n.samples)
data.cond.2 <- generate.random.subsamples(cond.2[, -1], n.samples)

dcloc.AB <- DCloc(data.cond.1, data.cond.2)
dcglob.AB <- DCglob(data.cond.1, data.cond.2)



perform.repetitions.local(combined,
                          n.repetitions,
                          n.samples,
                          path = opt$output.path,
                          file.name = "/n0_DCloc_")
perform.repetitions.global(combined,
                           n.repetitions,
                           n.samples,
                           opt$output.path,
                           "/n0_DCglob")

# local false discovery estimate ----------------
n0.expectation.local <-
  local.n0.estimation(n.repetitions,
                      opt$output.path,
                      "/n0_DCloc_",
                      d.threshold = opt$dthresh)
n0.local <- estimate.null.expectation(n0.expectation.local)
nAB.local <-
  local.nr.genes(as.data.frame(dcloc.AB), d.threshold = opt$dthresh)

fdr.local <-
  calculate.FDR(nAB.local, n0.local$mean, nrow(data.cond.1))

print(
  paste0(
    "False discovery rate for local topology analysis for d = ",
    opt$dthresh,
    " :",
    fdr.local
  )
)

# global false discovery estimate ----------------
n0.expectation.global <-
  global.n0.estimation(n.repetitions,
                       opt$output.path,
                       "/n0_DCglob",
                       p.threshold = opt$pthresh)
n0.global <- estimate.null.expectation(n0.expectation.global)
nAB.global <-
  global.nr.genes(as.data.frame(dcglob.AB), p.threshold = opt$dthresh)

fdr.global <-
  calculate.FDR(nAB.global, n0.global$mean, nrow(data.cond.1))

print(
  paste0(
    "False discovery rate for global topology analysis for p = ",
    opt$pthresh,
    " :",
    fdr.global
  )
)
