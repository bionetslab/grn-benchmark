# Use the official R base image
FROM r-base:latest

# Install required system libraries
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages('BiocManager', repos='https://cran.rstudio.com')"
RUN R -e "BiocManager::install(c('pcaMethods', 'multtest', 'GEOquery', 'affy', 'genefilter', 'GOstats', 'ath1121501.db'), force=TRUE)"
RUN R -e "install.packages(c('spatstat', 'igraph'), repos='https://cran.rstudio.com')"
RUN R -e "install.packages('DiffCorr', repos='https://cran.rstudio.com')"

# Copy the R script into the container
COPY script_arbitrary.R /usr/local/bin/script_arbitrary.R

# Set the working directory
WORKDIR /usr/local/bin

# Command to run the script with arguments
ENTRYPOINT ["Rscript", "script_arbitrary.R"]
