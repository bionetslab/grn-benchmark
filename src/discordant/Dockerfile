FROM rocker/r-ver:latest
RUN /rocker_scripts/install_tidyverse.sh

LABEL Maintainer="Karim Sakr <karim.m.sakr@fau.de>"

USER root

RUN apt-get update && apt-get install -y  gcc libglpk40

RUN R -e "install.packages(c('BiocManager', 'optparse', 'data.table'),dependencies=TRUE)"

RUN R -e "BiocManager::install(c('discordant', 'clusterProfiler', 'org.Hs.eg.db'))"

COPY . /discordant

WORKDIR /discordant