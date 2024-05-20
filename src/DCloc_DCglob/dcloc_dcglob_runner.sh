Rscript dcloc_dcglob.R -c data/out_CD8_exhausted.tsv -d data/out_Macrophages.tsv -o output/ --global --pthresh 0.1 --corrthresh 0.5
Rscript dcloc_dcglob.R -c data/out_CD8_exhausted.tsv -d data/out_Macrophages.tsv -o output/ --local --dthresh 0.25 --corrthresh 0.5

