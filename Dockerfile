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
		libssl-dev \
		zlib1g-dev \
		libcurl4-openssl-dev \ 
		openssl \
		default-jdk \
		lbzip2 \
		unzip \
		bzip2 \
		python \
		python2.7-dev \
		python-matplotlib \
		python-numpy \
       	python-pandas \
		tzdata \ 
		perl \
		wget \
		bcftools \
		curl

RUN apt-get install -y locales && locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

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
	conda install -c bioconda diamond && \
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
# Diamond Blast
#RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.10/diamond-linux64.tar.gz && \
#	tar xzf diamond-linux64.tar.gz
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

# Uniprot database
ADD https://github.com/iPlantCollaborativeOpenSource/docker-builds/releases/download/evolinc-I/uniprot_sprot.dmnd.gz /evolinc_docker/
RUN gzip -d /evolinc_docker/uniprot_sprot.dmnd.gz

# rFAM database
ADD https://de.cyverse.org/dl/d/12EF1A2F-B9FC-456D-8CD9-9F87197CACF2/rFAM_sequences.fasta /evolinc_docker/

# Biopython
RUN curl "https://bootstrap.pypa.io/pip/2.7/get-pip.py" -o "get-pip.py" && \
	curl "https://files.pythonhosted.org/packages/27/79/8a850fe3496446ff0d584327ae44e7500daf6764ca1a382d2d02789accf7/pip-20.3.4-py2.py3-none-any.whl" -o "pip-20.3.4-py2.py3-none-any.whl" && \
	curl "https://files.pythonhosted.org/packages/e1/b7/182161210a13158cd3ccc41ee19aadef54496b74f2817cc147006ec932b4/setuptools-44.1.1-py2.py3-none-any.whl" -o "setuptools-44.1.1-py2.py3-none-any.whl" && \
	curl "https://files.pythonhosted.org/packages/27/d6/003e593296a85fd6ed616ed962795b2f87709c3eee2bca4f6d0fe55c6d00/wheel-0.37.1-py2.py3-none-any.whl" -o "wheel-0.37.1-py2.py3-none-any.whl" && \
	python get-pip.py "pip-20.3.4-py2.py3-none-any.whl" "setuptools-44.1.1-py2.py3-none-any.whl" "wheel-0.37.1-py2.py3-none-any.whl" && \
	curl "https://files.pythonhosted.org/packages/ff/f4/0ce39bebcbb0ff619426f2bbe86e60bc549ace318c5a9113ae480ab2adc7/biopython-1.76.tar.gz" -o "biopython-1.76.tar.gz" && \
	tar xvf biopython-1.76.tar.gz && \
	cd biopython-1.76 && python setup.py build && python setup.py install && \
	cd ..

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
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
	tar -xvjf samtools-1.10.tar.bz2 
WORKDIR samtools-1.10 
RUN ./configure --prefix=/usr/bin && \
	make && \
	make install
WORKDIR /

#HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
	tar -xvjf htslib-1.17.tar.bz2 && \
	cd htslib-1.17 && \
	make && \
	cd ..

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
ENV PATH /htslib/:$PATH
ENV PATH /HAMR/hamr.py:$PATH
ENV PATH /HAMR/:$PATH

# Caching the sra data
RUN vdb-config --root -s /repository/user/cache-disabled="true"

# pamlinc wrapper script
ADD pamlinc.sh $BINPATH
RUN chmod +x $BINPATH/pamlinc.sh

ENTRYPOINT ["pamlinc.sh"]
