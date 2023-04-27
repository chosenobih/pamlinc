#!/bin/bash

#Add conda channels
conda config --add channels conda-forge 
conda config --add channels defaults  
conda config --add channels bioconda

#install conda packages
echo "installing packages required to run pamlinc"

conda install star==2.7.10a -c bioconda -y
conda install trimmomatic==0.35 -c bioconda -y
conda install gffread==0.12.1 -c bioconda -y
conda install subread==2.0.1 -c bioconda -y
conda install samtools==1.9 -c bioconda -y
conda install cufflinks==2.2.1 -c bioconda -y
conda install stringtie==2.1.5 -c bioconda -y
conda install picard==2.18.27 -c bioconda -y

echo "installation done"
