# PAMLINC: Parallel Annotation of Modified RNA and LincRNA

## Introduction

* PAMLINC is a workflow for processing raw RNA-Seq illumina data for rapid quantification of transcript abundance, RNA modifications and Long non-coding Intergenic ribonucleotides (lincRNA).
* PMALINC can process raw FASTQ files containing either paired-end or single-end reads. It can also process sequence read archive (SRA) from NCBI using a SRA ID. PAMLINC supports two reads aligner options: tophat2 and STAR but we recommend that users who are interested in using PAMLINC to annotate RNA modification select tophat2 as the aligner of choice due to alignment alogirthms compatibility.
* PAMLINC minimally requires the following input data:
  1. Reference genome (FASTA)
  2. Reference annotation (GTF/GFF3)
  3. RNA-Seq reads (FASTQ) - Paired end or Single end or NCBI SRA ID.

PAMLINC command line arguments and description
----------------------------------------------
| Argument      | Description                                                                                                 |
| ------------- |:-----------------------------------------------------------------------------------------------------------:|
| -g            | reference genome fasta file                                                                                 |
| -a            | reference genome annotation file                                                                            |
| -A            | reference genome annotation file type ("GTF" or "GFF3" supported. Include double quotation on command line) |
| -i            | index folder                                                                                                |
| -l            | library type  (library type can be fr-unstranded, fr-firststrand or fr-secondstrand)                        |
| -1            | read_1                                                                                                      |
| -2            | read_2                                                                                                      |
| -u            | single_reads                                                                                                |
| -o            | output directory                                                                                            |
| -S            | NCBI SRA-ID                                                                                                 |
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
2. [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html)

Running PAMLINC
-----------------------
Pull the container:  
```
singularity pull docker://
```  
