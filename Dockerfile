FROM ubuntu:18.04

LABEL This Dockerfile is for pamlinc. It is maintained by Chosen Obih <chosenobih@arizona.edu>
ENV DEBIAN_FRONTEND=noninteractive

USER root

RUN apt-get update && apt-get install -y g++ \
		build-essential \
	   	make \
		git \
		libcurl4 \
		libcurl4-openssl-dev \
		libssl-dev \
		libncurses5-dev \
		libsodium-dev \
		libmariadb-client-lgpl-dev \
		libbz2-dev \
		liblzma-dev \
		libssl-dev \
		zlib1g-dev \
		libcurl4-openssl-dev \ 
		openssl \
		default-jdk \
		lbzip2 \
		unzip \
		bzip2 \
		python2.7 \
		python2.7-dev \
		python-pip \
		python-matplotlib \
		python-numpy \
       	python-pandas \
		tzdata \ 
		perl \
		wget \
		bcftools \
		curl

RUN ldconfig
RUN apt-get install -y locales && locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

RUN pip2 install biopython==1.76

# Downlaod and install conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Conda channels
RUN conda config --add channels conda-forge && \
	conda config --add channels bioconda

# Conda packages
RUN conda install fastp==0.23.4 -c bioconda -y && \
	conda install biopython==1.76 -c anaconda -y && \
	conda install star==2.7.10a -c bioconda -y && \
	conda install trimmomatic==0.35 -c bioconda -y && \
	conda install gffread==0.12.1 -c bioconda -y && \
	conda install subread==2.0.1 -c bioconda -y && \
	conda install samtools==1.9 -c bioconda -y && \
	conda install cufflinks==2.2.1 -c bioconda -y && \
	conda install stringtie==2.1.5 -c bioconda -y && \
	conda install picard==2.18.27 -c bioconda -y && \
	conda install bowtie2==2.2.5 -c bioconda -y && \
	conda install tophat==2.1.1 -c bioconda -y && \
	conda install sra-tools==3.0.0 -c bioconda -y && \
	conda install numpy -y && \
	conda install pandas -y && \
	conda install last -c bioconda -y && \
	conda install diamond==0.9.10 -c bioconda && \
	conda install matplotlib-base -c conda-forge -y


# Required files
WORKDIR /

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && \
	unzip Trimmomatic-0.38.zip && rm Trimmomatic-0.38.zip && \
	cp -R /Trimmomatic-0.38/adapters /usr/bin/adapters
ENV ADAPTERPATH=/usr/bin/adapters

## evolinc-part-I
RUN git clone https://github.com/chosenobih/Evolinc-I.git
RUN cp -R /Evolinc-I /evolinc_docker
ENV BINPATH /usr/bin
WORKDIR /evolinc_docker

# Cufflinks
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Transdecoder
RUN wget -O- https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz | tar xzvf -

# Bedtools
RUN wget https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz && \
	tar xvf v2.25.0.tar.gz && \
	cd bedtools2-2.25.0 && make && \
	cd ..

# Bedops tool
RUN wget -O- https://github.com/bedops/bedops/releases/download/v2.4.16/bedops_linux_x86_64-v2.4.16.tar.bz2 | tar jxvf -
# cpan
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm URI/Escape.pm

# R libraries
RUN apt-get -y install ca-certificates software-properties-common gnupg2 gnupg1 gnupg && \
	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" && \
	apt-get install -y r-base && \
	Rscript -e 'install.packages("openssl", dependencies = TRUE,  repos="http://cran.rstudio.com/")' && \
	Rscript -e 'install.packages("splitstackshape", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
	Rscript -e 'install.packages("dplyr", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
	Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")};' && \
	Rscript -e 'BiocManager::install(c("Biostrings"));' && \
	Rscript -e 'install.packages("getopt", dependencies = TRUE, repos="http://cran.rstudio.com/");'

# Remove the existing symbolic link (if it exists)
RUN rm /usr/bin/python

# Create a symbolic link to make 'python' refer to 'python2.7'
RUN ln -sf /usr/bin/python2.7 /usr/bin/python

# Set Python 2.7 as the default python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1

# Uniprot database
ADD https://github.com/iPlantCollaborativeOpenSource/docker-builds/releases/download/evolinc-I/uniprot_sprot.dmnd.gz /evolinc_docker/
RUN gzip -d /evolinc_docker/uniprot_sprot.dmnd.gz && \
	chmod +r uniprot_sprot.dmnd

# rFAM database
ADD https://de.cyverse.org/dl/d/12EF1A2F-B9FC-456D-8CD9-9F87197CACF2/rFAM_sequences.fasta /evolinc_docker/

# CPC2
WORKDIR /evolinc_docker/CPC2-beta/libs/libsvm/
RUN tar xvf libsvm-3.22.tar.gz
WORKDIR libsvm-3.22
RUN make clean && make
WORKDIR /

# evolinc-part-I wrapper script
RUN chmod +x /evolinc_docker/evolinc-part-I.sh && cp /evolinc_docker/evolinc-part-I.sh $BINPATH

## HAMR
RUN git clone https://github.com/chosenobih/HAMR.git && \
	chmod +x /HAMR/hamr.py && cp /HAMR/hamr.py $BINPATH && \
	cp -R /HAMR/models /usr/bin/hamr_models
ENV HAMR_MODELS_PATH=/usr/bin/hamr_models

# samtools & 
#RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
#	tar -xvjf samtools-1.10.tar.bz2 && \
# 	cd samtools-1.10 && \
#	./configure --prefix=/usr/bin && \
#	make && \
#	make install
#WORKDIR /

#HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
	tar -xvjf htslib-1.17.tar.bz2 && \
	cd htslib-1.17 && \
	./configure --prefix=/usr/bin && \
	make && \
	make install
WORKDIR /

#GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip && \
	unzip gatk-4.2.2.0.zip && rm gatk-4.2.2.0.zip
ENV PATH="/gatk-4.2.2.0:${PATH}"

# Setting paths to all the softwares
ENV PATH /evolinc_docker/cufflinks-2.2.1.Linux_x86_64/:$PATH
ENV PATH /evolinc_docker/TransDecoder-2.0.1/:$PATH
ENV PATH /evolinc_docker/ncbi-blast-2.4.0+/bin/:$PATH
ENV PATH /evolinc_docker/bedtools2-2.25.0/bin/:$PATH
ENV PATH /evolinc_docker/bin/:$PATH
ENV PATH /evolinc_docker/CPC2-beta/bin/:$PATH
ENV PATH /evolinc_docker/:$PATH
ENV PATH /usr/bin/:$PATH
ENV PATH /HAMR/hamr.py:$PATH
ENV PATH /HAMR/:$PATH

# Caching the sra data
RUN vdb-config --root -s /repository/user/cache-disabled="true"

# pamlinc wrapper script
ADD pamlinc.sh $BINPATH
RUN chmod +x $BINPATH/pamlinc.sh

ENTRYPOINT ["pamlinc.sh"]