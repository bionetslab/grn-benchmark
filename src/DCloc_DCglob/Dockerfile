
FROM rocker/r-ver:4.2.3


RUN apt-get update && apt-get install -y \
    bash \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /home/analysis

RUN R -e "install.packages(c('BiocManager', 'optparse', 'data.table', 'gplots', 'tidyverse', 'igraph'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('affy')"
RUN R -e "install.packages('tidyverse')"

COPY . /home/analysis

CMD cd /home/analysis


