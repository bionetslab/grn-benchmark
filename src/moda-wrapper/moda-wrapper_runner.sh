#!/usr/bin/env bash

Rscript ./R/hierarchical-clustering.R ./R/dataset/Macrophages/ ./R/dataset/TCells/ ./R/output/ -i1 Macrophages -i2 TCells -c 2 -s 0.1 -t 0.1
