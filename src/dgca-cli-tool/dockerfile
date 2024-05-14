# Use an R base image with version 4.0.5 from Rocker
FROM bioconductor/bioconductor_docker:RELEASE_3_17

# Install necessary Linux dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev && \
    apt-get clean

# Copy the dependencies installation script into the container
COPY install_dependencies.R /usr/local/src/scripts/install_dependencies.R

# Run the install dependencies script
RUN Rscript /usr/local/src/scripts/install_dependencies.R

# Copy the R scripts into the container
COPY dgca.R /usr/local/src/scripts/dgca.R
COPY enrichment_map.R /usr/local/src/scripts/enrichment_map.R

# Make the R scripts executable
RUN chmod +x /usr/local/src/scripts/dgca.R
RUN chmod +x /usr/local/src/scripts/enrichment_map.R

# Set the working directory to the scripts directory
WORKDIR /usr/local/src/scripts

# Set entrypoint to allow general Rscript execution
ENTRYPOINT ["Rscript"]

