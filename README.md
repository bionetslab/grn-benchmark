# BIONETS Hackathon repository

This is the Github project page for the BIONETS project / hackathon in the summer term 2024.

This wiki may eventually be made public.

# How to contribute a project
## Preparation
1. Pick next available project and assign it to yourself (see supplied list of references)
2. Read paper and update summary table with missing information
3. Find the code or software online
4. Check if tutorials exist by the authors and update summary sheet.
5. Check if tutorials exist by external authors
6. Find which data has been used by the authors and how it was preprocessed.
7. Update the data set information sheet with the relevant information.
8. If available, find the code used for the preprocessing of the data.
9. Find the settings the authors used in the publication to generate the figures.
10. Create a new subfolder using the tool name (all lower case, hyphenated) and implement the tool

## Implementation:
1. Install software and create log file according to this file: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-09406-4/MediaObjects/41467_2019_9406_MOESM1_ESM.pdf
2. Check if Minimal working example exists and run if available
3. Check if code exists for creating the figures in the article
4. Attempt to replicate the examples/figures shown in the if not available using the data set(s) supplied in the study.
5. Create script allowing the execution of the tool using all reference data sets.
7. Document method parameters, inputs and outputs.
8. Create Docker container
9. Supply yaml file for a conda environment.
10. Push code to github repo
11. Create markdown README in project folder
12. Wiki page with methodology, rationale, parameters, etc.
13. Update secondary evaluation criteria list

## A note on the scripts.
1. Follow the supplied specifications regarding parameters, output folder strucutre, etc. 
2. A script should allow you to execute a tool using one dataset with one parameter setting
3. If you want to test multiple parameter settings please create a wrapper script which calls the script with the relevant parameters.
4. If the tool is a commandline tool itself, it is not necessary to wrap the tool again.
5. If the tool/library is written in R, the script should be callable using R: ```Rscript tool-name.R -p p1 -q p2```.


## Troubleshooting:
1. If the installation fails, troubleshoot issues, double check if someone else is able to install it on their computer. (Especially with R, you sometimes need to install new system libraries, therefore you can use conda or docker, ...)
2. If the execution of the software fails, troubleshoot the issue (memory error) and try to fix it. Otherwise report issue.
3. If the running of the tool takes unreasonable time (e.g. >2h for a small example dataset) try running over night and report the run time.
4. If other problems occur, please document your issues. 


## README
For every tool there should be a README with
1. Brief description of the tool
2. Reference to the publication
3. Installation instructions, or relevant links to the instructions if there were no issues you encountered.
4. Copy-and-pastable execution instructions using example data.
5. Explanation of the relevant parameters
6. Input file format specification
7. Output file format specification
8. Explanation and interpretation of the output
9. Recommended hyperparameters by the authors
10. Hyperparameter recommendations for optimization (more instructions will follow)
11. Other necessary information

In general, the more difficult it was to install, execute or interpret the results of the tool, the more information needs to be supplied in the README.md file.
