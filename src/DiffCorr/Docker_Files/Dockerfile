# Use an R base image as the base image
FROM rocker/r-ver:4.1.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libgmp-dev \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN install2.r --error \
    BiocManager

# Install R packages from CRAN and Bioconductor
RUN R -e "BiocManager::install(c('affy', 'genefilter', 'igraph', 'spatstat', 'cluster', 'DiffCorr'))"


# Copy your R scripts and any other necessary files into the container
COPY downstream_macro.R downstream_macro.R
COPY script_macro_cd8.R script_macro_cd8.R
COPY downstream_gse.R downstream_gse.R
COPY script_gse.R script_gse.R
COPY all_aml.R all_aml.R

# Set the working directory
WORKDIR /app

# Run your R script when the container starts
CMD ["Rscript", "script_macro_cd8.R"]