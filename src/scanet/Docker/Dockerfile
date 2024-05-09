# Use rocker/r-ver:latest as the base image
FROM rocker/r-ver:latest

# Install Python 3.11.9 and pip
RUN apt-get update && \
    apt-get install -y python3.11 python3.11-distutils python3-pip && \
    ln -sf /usr/bin/python3.11 /usr/bin/python3 && \
    ln -sf /usr/bin/python3.11 /usr/bin/python && \
    apt-get install -y python3-tk
    
RUN apt-get install -y libglpk40


# Install system dependencies for Python packages
RUN apt-get install -y \
    python3.11-dev \
    libicu-dev \
    libpcre2-dev \
    libdeflate-dev \
    liblzma-dev \
    libbz2-dev

# Set the working directory in the container to /app
WORKDIR /app


# Update pip and install required Python packages
RUN pip install --upgrade pip && \
    pip install scanet numpy==1.23.4 rpy2 pyvis drugstone tk threadpoolctl==3.1.0

# Install R packages
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('BioNERO', 'R.utils', 'testthat', 'ddpcr', 'tibble', 'DESeq2'))"    
# Optionally, copy your application code into the container
# COPY . /app

# Copy only the necessary files into the container
# It's recommended to use a .dockerignore file to exclude files not needed
COPY . /app

RUN rm -rf /usr/local/lib/python3.11/dist-packages/scanet/grn.py
RUN cp /app/grn.py /usr/local/lib/python3.11/dist-packages/scanet/grn.py

# Make port 80 available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME World

# Instead of specifying a default command, set up the container to run any script as needed
CMD ["bash"]
