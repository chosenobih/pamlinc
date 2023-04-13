# PAMLINC: Parallel Annotation of Modified RNA and LincRNA

## Introduction

* PAMLINC is a workflow for processing raw RNA-Seq illumina data for rapid quantification of transcript abundance, RNA modifications and Long non-coding Intergenic ribonucleotides (lincRNA).
* PMALINC can process raw FASTQ files containing either paired-end or single-end reads. It can also process sequence read archive (SRA) from NCBI using a SRA ID. PAMLINC supports two reads aligner options: tophat2 and STAR but we recommend that users who are interested in using PAMLINC to annotate RNA modification select tophat2 as the aligner of choice due to alignment alogirthms compatibility.
* PAMLINC minimally requires the following input data:
  1. Reference genome (FASTA)
  2. Reference annotation (GTF/GFF3)
  3. RNA-Seq reads (FASTQ) - Paired end or Single end or NCBI SRA ID.
* Optional files:
    The -i flag allows users to provide a reference genome index folder which should contain the genome index files for either bowtie2, STAR or both,           depending on the user's aligner of choice. The STAR index folder should be named 'star_index' and it should be a subdirectory of the reference genome       index folder provided. PAMLINC automatically generates the reference genome index files for both bowtie2 and STAR when it is not provided by the user       but this increases the run time of PAMLINC.


PAMLINC command line arguments and description
------------------------------------------------------------------------------------------------------------------------------
| Argument      | Description                                                                                                 |
| ------------- |:-----------------------------------------------------------------------------------------------------------:|
| -g            | reference genome fasta file                                                                                 |
| -a            | reference genome annotation file                                                                            |
| -i            | reference genome index folder                                                                               |
| -l            | library type  (library type can be fr-unstranded, fr-firststrand or fr-secondstrand)                        |
| -1            | read_1                                                                                                      |
| -2            | read_2                                                                                                      |
| -u            | single_reads                                                                                                |
| -o            | output directory                                                                                            |
| -S            | NCBI SRA-ID #upper case S                                                                                   |
| -p            | number of threads                                                                                           |
| -q            | activate transcript abundance quantification option                                                         |
| -t            | selects tophat2 as aligner of choice                                                                        |
| -s            | selects STAR as aligner of choice                                                                           |
| -y            | type of read (single end or paired end) #denoted as "SE" or "PE", include double quotation on command line  |
| -b            | reads_mismatches (% reads mismatches to allow)                                                              |
| -m            | activates RNA modification annotation option                                                                |
| -e            | activate lincRNA annotation option                                                                          |
| -k            | feature_type #Feature type (Default is exon)                                                                |
| -r            | gene attribute (Default is gene_id)                                                                         |
| -n            | strandedness (Default is 0 (unstranded), 1 (stranded), 2 (reversely stranded)                               |

Required dependecies
--------------------
1. Linux-based computer, server or cluster.
2. [Singularity v3.8 and above](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html)
3. [Conda](https://conda.io/projects/conda/en/stable/user-guide/install/download.html)

Running PAMLINC
-----------------------
Download and install singularity and conda

```
#clone the repo:  
https://github.com/chosenobih/pamlinc.git
```  
```
#run setup_script.sh
bash setup_script.sh
```
```
#download genome file from CyVerse data store
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sbicolor.fa
```
```
#download genome annotation file from CyVerse data store
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sbicolor.gff3
```

Running pamlinc in paired-end mode
```
#download sample data from CyVerse data store:
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sample_1_R1.fastq.gz
wget https://data.cyverse.org/dav-anon/iplant/home/chosen/pamlinc_files/sample_1_R2.fastq.gz
```
```
#run pamlinc to annotate RNA modification, identify lincRNA and quantify transcript abundance with paired fastq.gz files. These options can be turned on or off using different flags.
bash pamlinc_main.sh -a sbicolor.gff3 -g sbicolor.fa -o pamlinc_result -y "PE" -p 6 -l fr-secondstrand -q -e -m -k exon -r gene_id -n 0 -d 12 -t -1 sample_1_R1.fastq.gz -2 sample_1_R2.fastq.gz
```
