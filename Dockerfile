FROM ubuntu:18.04

LABEL This Dockerfile is for pamlinc. It is maintained by Chosen Obih <chosenobih@arizona.edu>
ENV DEBIAN_FRONTEND=noninteractive

USER root

RUN apt-get update && apt-get install -y \
	    	make \
		git \
		zlib1g-dev \
		openssl \
		lbzip2 \
		bzip2 \
		python \
		python2.7-dev \
		perl \
		wget \
		curl

#RUN add-apt-repository -y ppa:openjdk-r/ppa
#RUN apt-get update && apt-get install -y openjdk-8-jdk libbz2-dev liblzma-dev

# Downlaod and install conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Conda channels
RUN conda config --add channels conda-forge 
RUN conda config --add channels defaults  
RUN conda config --add channels bioconda
RUN conda config --add channels r

# Conda packages
RUN conda install fastp==0.23.4 -c bioconda -y
RUN conda install star==2.7.10a -c bioconda -y
RUN conda install trimmomatic==0.35 -c bioconda -y
RUN conda install gffread==0.12.1 -c bioconda -y
RUN conda install subread==2.0.1 -c bioconda -y
RUN conda install samtools==1.9 -c bioconda -y
RUN conda install cufflinks==2.2.1 -c bioconda -y
RUN conda install stringtie==2.1.5 -c bioconda -y
RUN conda install picard==2.18.27 -c bioconda -y

# Required files
WORKDIR /
RUN wget https://raw.githubusercontent.com/chosenobih/pamlinc/main/TruSeq3-PE.fa
RUN wget https://raw.githubusercontent.com/chosenobih/pamlinc/main/TruSeq3-SE.fa
RUN wget https://github.com/chosenobih/pamlinc/blob/main/GenomeAnalysisTK.jar
ENV PATH /TruSeq3-PE.fa/:$PATH
ENV PATH /TruSeq3-SE.fa/:$PATH
ENV PATH /GenomeAnalysisTK.jar/:$PATH

ENV BINPATH /usr/bin
ENV LC_ALL C 

# Wrapper script
ADD pamlinc_main.sh $BINPATH
RUN chmod +x $BINPATH/pamlinc_main.sh

ENTRYPOINT ["pamlinc_main.sh"]
