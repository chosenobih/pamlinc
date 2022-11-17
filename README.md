# PAMLINC: High-throughput Distributed Computing Pipeline for Processing and Analyzing RNA-Seq Data

## Introdcution

* PAMLINC is a workflow for processing raw RNA-Seq illumina data for rapid quantification of transcript abundance, RNA modifications and Long non-coding Intergenic ribonucleotides (lincRNA).
* PMALINC can process raw FASTQ files containing either paired-end or single-end reads. It can also process seuqnece read archive (SRA) from NCBI using a SRA ID.
* PAMLINC supports two reads aligner options: tophat2 and STAR but we recommend that users who are interested in using PAMLINC to annotate RNA modification annotation select tophat2 as the aligner of choice.
* PAMLINC minimally requires the following input data:
  1. Reference genome (FASTA)
  2. Reference annotation (GTF/GFF3)
  3. RNA-Seq reads (FASTQ) - Paired end or Single end or NCBI SRA ID.

Pipeline programs and version
-----------------------------
| Command Line  | Description    |
  Argument
| ------------- |:-------------:|
|               |           |
| TopHat2       |          |
| Samtools      |           |
| Picard        |         |
| Gakt          |          |
| HAMR          |            |
| Evolinc-i     |          |
| Salmon        |          |
| subread       |          |
| sra-tools     |            |

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
