1. Installation
2. System requirements


 Installation

 NOTE that this version of Diniro has been generated for local usage.
 The webtool is in a docker image that is generated using docker-compose. 
 
You need first to download the docker for desktop following this link:
https://docs.docker.com/compose/

 Once this is done, all you need to do is:
 
1. Using cmd clone the repository from github

2. After you clone it, open the cloned repo directly and open cmd from there

3. Make sure the the docker desktop is open and build your docker through: docker compose build 

4. If you didn't face any problem run: docker-compose up

5.  Go to: http://0.0.0.0:8025/

6. If you faced any errors in build the docker make sure to check the requirements.txt file and the dependancies 

You can then use Diniro LOCALY.
 
Note: 
To run the tool correctly, you need to download these files
1- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather" from this link:  https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/
2- "motifs-v9-nr.hgnc-m0.001-o0.0.tbl" from this link: https://resources.aertslab.org/cistarget/motif2tf/
and add it to this directory: 
diniro\Docker\test_interface\static\GRN_files\human

I couldn't upload it GitHub Repo.
 
 
 Documentation

 See https://exbio.wzw.tum.de/diniro// documentation or install the program.
 
 
 Hardware requirements:
 The program requires a standard computer with enough RAM to support the in-memory operations.
 

 OS requirements:
 The program has only been tested on:
 Linux: Ubuntu: 16.04
 Windows 11
 










