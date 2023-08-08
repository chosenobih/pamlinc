# PAMLINC: Parallel Annotation of Modified RNA and LincRNA

![pamlinc_workflow drawio](https://github.com/chosenobih/pamlinc/assets/50637858/2465bfb4-2c3f-4f70-a69d-041288caf2df)

## Introduction

* PAMLINC is a workflow for processing raw RNA-Seq illumina data for rapid quantification of transcript abundance, annotation of RNA modifications and Long non-coding Intergenic ribonucleotides (lincRNA).
* PAMLINC can process raw FASTQ files containing either paired-end or single-end reads. It can also process sequence read archive (SRA) from NCBI using a SRA ID. PAMLINC supports two reads trimming options: trimmomatic and fastp. PAMLINC also supports two reads aligner options: TopHat2 and STAR, but users who are interested in using PAMLINC to annotate RNA modification must select TopHat2 as the aligner of choice. This is due to alignment algorithms compatibility.
* PAMLINC minimally requires the following input data:
  1. Reference genome (FASTA)
  2. Reference annotation (GTF/GFF3)
  3. RNA-Seq reads (FASTQ) - Paired end or Single end or NCBI SRA ID.
* Optional files:
    The -i flag allows users to provide a reference genome index folder which should contain the genome index files for either bowtie2, STAR or both, depending on the user's aligner of choice. The STAR index folder should be named 'star_index' and it should be a subdirectory of the reference genome index folder provided. PAMLINC automatically generates the reference genome index files for both bowtie2 and STAR when it is not provided by the user, but this increases the run time of PAMLINC.
* Output: When run with all three flags (-m, -e, -q), PAMLINC generates four folders in the output directory. The folders are <sample_name>_HAMR, <sample_name>_lincRNA, intermediate_files, and transcript_abund_quant


PAMLINC command line arguments and description
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
| Argument      | Description                                                                                                                                                                       |
| ------------- |:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| -g            | reference genome fasta file                                                                                                                                                       |
| -a            | reference genome annotation file                                                                                                                                                  |
| -i            | reference genome index folder                                                                                                                                                     |
| -l            | library type  (library type can be fr-unstranded, fr-firststrand or fr-secondstrand). If you are not sure of the library type of your reads, you can infer it using salmon.       |
| -1            | read_1                                                                                                                                                                            |
| -2            | read_2                                                                                                                                                                            |
| -u            | single_read                                                                                                                                                                       |
| -o            | output directory                                                                                                                                                                  |
| -S            | NCBI SRA-ID #upper case S                                                                                                                                                         |
| -p            | number of threads                                                                                                                                                                 |
| -q            | activate transcript abundance quantification option                                                                                                                               |
| -t            | selects TopHat2 as aligner of choice                                                                                                                                              |
| -s            | selects STAR as aligner of choice                                                                                                                                                 |
| -w            | use trimmomatic for trimming input reads                                                                                                                                          |
| -x            | use fastp for trimming input reads                                                                                                                                                |
| -c            | sjdbOverhang #Required for STAR aligner. Value is dependent on the read length of your fastq files (read length minus 1)                                                          |
| -f            | genomeSAindexNbases #Required for STAR aligner. Default value is 14 but it should be scaled downed for small genomes, with a typical value of min(14, log2(GenomeLength)/2 - 1)   |
| -y            | type of read (single end or paired end) #Denoted as "SE" or "PE", include double quotation on command line                                                                        |
| -d            | reads_mismatches (% reads mismatches to allow) #Required for TopHat2 aligner                                                                                                      |
| -m            | activates RNA modification annotation option                                                                                                                                      |
| -e            | activate lincRNA annotation option                                                                                                                                                |
| -E            | run evolinc_i with mandatory files or both mandatory and optional files. Denoted as "M" or "MO", include double quotation on the command line                                     |
| -T            | </path/to/transposable Elements file>                                                                                                                                             |
| -C            | </path/to/CAGE RNA file>                                                                                                                                                          |
| -D            | </path/to/known lincRNA file>                                                                                                                                                     |
| -k            | feature_type #Feature type (Default is exon) #Required for featurecount                                                                                                           |
| -r            | gene attribute (Default is gene_id) #Required for featurecount                                                                                                                    |
| -n            | strandedness (Default is 0 (unstranded), 1 (stranded), 2 (reversely stranded) #Required for featurecount                                                                          |


Required dependecies
--------------------
1. Linux-based computer, server, or cluster.
2. [Docker](https://docs.docker.com/engine/install/)

Running PAMLINC
-----------------------

```
#pull docker image for PAMLINC:  
docker pull chosenobih/pamlinc:v0.1
```  

```
#download genome file for sorghum bicolor from CyVerse data store
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sbicolor.fa
```
```
#download genome annotation file for sorghum bicolor from CyVerse data store
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sbicolor.gff3
```
```
#download sbicolor genome index files for bowtie-2 from google drive. Copy the link and paste it into your browser.
https://drive.google.com/drive/folders/1EazCZ1__K7DCKbOOYoQ8jeqj3Tpb9KrR?usp=sharing
```
Running pamlinc in paired-end mode with a fastq file
```
#download sample data from CyVerse data store:
#R1
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sample_1_R1.fastq.gz
```
```
#download sample data from CyVerse data store:
#R2
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sample_1_R2.fastq.gz
```

```
#run pamlinc to annotate RNA modification, identify lincRNA and quantify transcript abundance with paired fastq.gz files. These options can be turned on or off using different flags.
docker run --rm -v $(pwd):/working-dir -w /working-dir chosenobih/pamlinc:v0.1 -a sbicolor.gff3 -g sbicolor.fa -o pamlinc_result_PE -y "PE" -p 6 -w -l fr-secondstrand -q -e -E "M" -m -k exon -r gene_id -n 0 -d 12 -t -1 sample_1_R1.fastq.gz -2 sample_1_R2.fastq.gz -i index_folder
```
```
#download paired-end sample run output from google drive. Copy the link and paste it into your browser.
https://drive.google.com/drive/folders/1A4CuPfrvX0oBw2cmqn8wOBI9_rMWdds0?usp=sharing
```

Running pamlinc in single-end mode with a fastq file
```
#download sample data from CyVerse data store:
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sample_1_SE.fastq.gz
```
```
#run pamlinc to annotate RNA modification, identify lincRNA and quantify transcript abundance with paired fastq.gz files. These options can be turned on or off using different flags.
docker run --rm -v $(pwd):/working-dir -w /working-dir chosenobih/pamlinc:v0.1 -a sbicolor.gff3 -g sbicolor.fa -o pamlinc_result_SE -y "SE" -p 6 -w -l fr-secondstrand -q -e -E "M" -m -k exon -r gene_id -n 0 -d 12 -t -u sample_1_SE.fastq.gz -i index_folder
```
```
#download single-end sample run output from google drive. Copy the link and paste it into your browser.
https://drive.google.com/drive/folders/1iLB7Gec9qB6Sv2TyFlHM_GG3HQ7wPWa9?usp=sharing
```

Running pamlinc in paired-end mode with an SRA-ID
```
#run pamlinc to annotate RNA modification, identify lincRNA and quantify transcript abundance with paired fastq.gz files. These options can be turned on or off using different flags.
docker run --rm -v $(pwd):/working-dir -w /working-dir chosenobih/pamlinc:v0.1 -a sbicolor.gff3 -g sbicolor.fa -o pamlinc_result_SRA-ID_PE -t -y "PE" -p 6 -w -S SRR18095197 -l fr-secondstrand -q -e -E "M" -m -k exon -r gene_id -n 0 -d 12 -i index_folder
```

Running pamlinc in single-end mode with an SRA-ID
```
#run pamlinc to annotate RNA modification, identify lincRNA and quantify transcript abundance with paired fastq.gz files. These options can be turned on or off using different flags.
docker run --rm -v $(pwd):/working-dir -w /working-dir chosenobih/pamlinc:v0.1 -a sbicolor.gff3 -g sbicolor.fa -o pamlinc_result_SRA-ID_SE -t -y "SE" -p 6 -w -S SRR10376271 -l fr-secondstrand -q -e -E "M" -m -k exon -r gene_id -n 0 -d 12 -i index_folder
```
