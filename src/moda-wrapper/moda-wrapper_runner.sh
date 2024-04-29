#!/usr/bin/env bash

Rscript ./R/hierarchical-clustering.R ./dataset/Macrophages/ ./dataset/TCells/ output/ -i1 Macrophages -i2 TCells -c 2 -s 0.1 -t 0.1
