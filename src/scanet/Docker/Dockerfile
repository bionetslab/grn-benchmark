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


WORKDIR /app


# Install Python packages
RUN pip install --upgrade pip && \
    pip install scanet numpy==1.23.4 rpy2 pyvis drugstone tk threadpoolctl==3.1.0

# Install R packages
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('BioNERO', 'R.utils', 'testthat', 'ddpcr', 'tibble', 'DESeq2'))"    

COPY . /app

RUN rm -rf /usr/local/lib/python3.11/dist-packages/scanet/grn.py
RUN cp /app/grn.py /usr/local/lib/python3.11/dist-packages/scanet/grn.py

EXPOSE 80

ENV NAME World

CMD ["bash"]
