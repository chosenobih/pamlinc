#!/bin/bash

#create new conda environment
conda create -n pamlinc_env python=3.10.6
conda activate pamlinc_env

#Add conda channels
conda config --add channels conda-forge 
conda config --add channels defaults  
conda config --add channels bioconda

#install conda packages
conda install star==2.7.10a -c bioconda
conda install trimmomatic==0.35 -c bioconda
conda install gffread==0.12.1 -c bioconda
conda install subread==2.0.1 -c bioconda
conda install samtools==1.9 -c bioconda
conda install cufflinks==2.2.1 -c bioconda
conda install stringtie==2.1.5 -c bioconda
conda install picard==2.18.27 -c bioconda
