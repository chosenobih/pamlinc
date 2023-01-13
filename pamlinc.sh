#!/bin/bash
#Chosen Obih
#Script to process and analyze RNA-Seq data for transcript abundance, modified RNA and lincRNA

usage() {
      echo ""
      echo "Usage : sh $0 -g <reference_genome>  -a <reference_annotation> -A <reference_transcriptome> -i <index_folder> -l lib_type {-1 <left_reads> -2 <right_reads> | -u <single_reads> | -S <sra_id>} -o <output_folder for pipeline files> -p num_threads -d reads_mismatches -t tophat -s star -c sjdbOverhang -f genomeSAindexNbases -q transcript_abundance_quantification -e evolinc_i -m HAMR -r gene_attribute -n strandedness -k feature_type"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########
  -g <reference genome fasta file>
  -a <reference genome annotation>
  -A <reference_transcriptome>
  -i <index_folder>
  -l library type #note that this is a lower case L
  -1 <reads_1>
               # Ends with R1 and is in the same order as reverse reads
  -2 <reads_2>
               # Ends with R2, must be present, and is in the same order as forward reads
  -u <single_reads> # Do not use single reads along with paired end reads
  -o </path/to/ pipeline output folder>
  -S SRA ID # one SRA ID or multiple SRA ID's in a .txt file with one 
  -p number of threads
  -q transcript abundance quantification
  -t tophat2 mapping #needed if you want to run HAMR
  -s star mapping #deactivates tophat2 and HAMR
  -c sjdbOverhang #value is dependent on the read length of your fastq files (read length minus 1)
  -f genomeSAindexNbases #default value is 14 but it should be scaled downed for small genomes, with a typical value of min(14, log2(GenomeLength)/2 - 1)
  -y type of reads (single end or paired end) #denoted as "SE" or "PE", include double quotation on command line
  -d reads_mismatches (% reads mismatches to allow. Needed for tophat2)
  -m HAMR
  -e evolinc_i
  -k feature type #Feature type (Default is exon)
  -r gene attribute (Default is gene_id)
  -n strandedness (Default is 0 (unstranded), 1 (stranded), 2 (reversely stranded)
EOF
    exit 0
}

star=0
tophat=0
referencegenome=0
referenceannotation=0
transcriptome=0

while getopts ":g:a:A:i:l:1:2:u:o:S:p:d:k:r:n:htsqemy:" opt; do
  case $opt in
    g)
    referencegenome=$OPTARG # Reference genome file
     ;;
    a)
    referenceannotation=$OPTARG # Reference genome annotation
     ;;
    A)
    transcriptome=$OPTARG # Reference genome transcriptome (transcript.fa)
     ;;
    i)
    index_folder=$OPTARG # Input folder
     ;;
    l)
    lib_type=$OPTARG # Library type (lib-type can be fr-unstranded, fr-firststrand or fr-secondstrand. If you are not sure of the library type of your reads, you can infer it using salmon.)
     ;;
    1)
    left_reads+=("$OPTARG") # Left reads
     ;;
    2)
    right_reads=("$OPTARG") # Right reads
     ;;
    u)
    single_reads+=("$OPTARG") # single end reads
     ;;
    o)
    pipeline_output=$OPTARG # pipeline output files
     ;;
    S)
    sra_id=$OPTARG # SRA ID or SRA ID's in a file
     ;;
    p)
    num_threads=$OPTARG # Number of threads
     ;;
    d)
    reads_mismatches=$OPTARG # Number of mismatches to allow in tophat2 run
     ;;
    q)
    transcript_abun_quant=$OPTARG # transcript abundance quantification
     ;;
    m)
    HAMR=$OPTARG # HAMR
     ;;
    e)
    evolinc_i=$OPTARG # evolinc_i
     ;;
    s)
    star=$OPTARG # star
     ;;
    c)
    sjdbOverhang=$OPTARG # need for star genome index building. Value is dependent on the read length of your fastq files (read length minus 1)
     ;;
    f)
    genomeSAindexNbases=$OPTARG # needed for star genome index building. Default value is 14 but it should be scaled downed for small genomes, with a typical value of min(14, log2(GenomeLength)/2 - 1)
     ;;
    t)
    tophat=$OPTARG # tophat
     ;;
    k) 
    feature_type=$OPTARG # Feature type (Default is exon)
     ;;
    r) 
    gene_attribute=$OPTARG # (Default is gene_id)
     ;;
    n)
    strandedness=$OPTARG # (Default is 0 (unstranded), 1 (stranded), 2 (reversely stranded))
     ;;
    y)
    seq_type=$OPTARG # Type of Sequence data (SE or PE. Mainly needed for SRA and featurecounts)
     ;;
    h)
    usage
     exit 1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Create the output directory
if [ ! -d "$pipeline_output" ]; then
  mkdir $pipeline_output
  elif [ -d "$pipeline_ouput" ]; then
  rm -r $pipeline_output; mkdir $pipeline_output
fi

###################################################################################################################
# # Check reference genome annotation file type and convert to .gtf if .gff file was supplied by user.
###################################################################################################################
if (grep -q -E 'transcript_id | gene_id' $referenceannotation); then
    echo "$referenceannotation is in .gtf format"
    else
    gffread $referenceannotation -T -o ref_annotation.gtf
    referenceannotation=ref_annotation.gtf
fi

paired_fastq_gz()
{
    filename=$(basename "$f" ".fastq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')

          if [ "$seq_type" == "PE" ]; then
          echo "######################"
          echo "Trimming input read(s)"
          echo "######################"
          
          echo "trimmomatic PE -threads $num_threads ${filename}.fastq.gz ${filename2}.fastq.gz ${filename3}_1P.fastq.gz ${filename3}_1U.fastq.gz ${filename3}_2P.fastq.gz ${filename3}_2U.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
          trimmomatic PE -threads $num_threads ${filename}.fastq.gz ${filename2}.fastq.gz ${filename3}_1P.fastq.gz ${filename3}_1U.fastq.gz ${filename3}_2P.fastq.gz ${filename3}_2U.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2 LEADING:3 TRAILING:3 MINLEN:36
          fi
          
          if [ "$tophat" != 0 ] && [ "$star" == 0 ]; then
          echo "###################################"
          echo "Running tophat2 in paired end mode"
          echo "###################################"
          
          if [ "$seq_type" == "PE" ]; then
          echo "tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_fwd_tophat -G $referenceannotation $fbname ${filename3}_1P.fastq.gz,${filename3}_1U.fastq.gz"
          singularity run --cleanenv tophat_2.1.1--py27_3.sif tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_fwd_tophat -G $referenceannotation $fbname ${filename3}_1P.fastq.gz,${filename3}_1U.fastq.gz
          echo "tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_rev_tophat -G $referenceannotation $fbname ${filename3}_2P.fastq.gz,${filename3}_2U.fastq.gz"
          singularity run --cleanenv tophat_2.1.1--py27_3.sif tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_rev_tophat -G $referenceannotation $fbname ${filename3}_2P.fastq.gz,${filename3}_2U.fastq.gz
          fi
          fi

}

#############################################################################################################################################################################################################################
# # Check if user supplied index folder for salmon, bowtie2 and star indexes. Build the indexes if not supplied.
#############################################################################################################################################################################################################################

if [ "$tophat" != 0 ] && [ "$star" == 0 ]; then
  if [ ! -z "$index_folder" ]; then
    for i in $index_folder/*; do
      mv -f $index_folder/*bt2l .
      fbname=$(basename "$i" .bt2l | cut -d. -f1)
    done
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
    echo "##########################################"
    echo "Building reference genome index for Tophat"
    echo "##########################################"
    echo "bowtie2-build --threads -f $referencegenome ref_genome"
    bowtie2-build -f $referencegenome ref_genome
    echo "fbname=$(basename "ref_genome" .bt2l | cut -d. -f1)"
    fbname=$(basename "ref_genome" .bt2l | cut -d. -f1)
  fi
elif [ "$tophat" == 0 ] && [ "$star" != 0 ]; then
  if [ ! -z "$index_folder" ]; then
      mv -f $index_folder/star_index/ .
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
    echo "########################################"
    echo "Building reference genome index for STAR"
    echo "########################################"
    echo "STAR --runThreadN $num_threads --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles $referencegenome --sjdbGTFfile $referenceannotation --sjdbOverhang 100 --genomeSAindexNbases 13"
    STAR --runThreadN $num_threads --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles $referencegenome --sjdbGTFfile $referenceannotation --sjdbOverhang 100 --genomeSAindexNbases 13
  fi
fi

if [ "$transcript_abun_quant" != 0 ]; then
  if [ ! -z "$index_folder" ]; then
      mv -f $index_folder/salmon_index/ .
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
  echo "####################################"
  echo "Building reference genome for salmon"
  echo "####################################"
  echo "salmon index -i salmon_index -t $transcriptome"
  salmon index -t $transcriptome -i salmon_index
  fi
fi
# ############################################################################################################################################################################################################################
# # Check that the input fastq files has the appropriate extension and then trim reads, align the reads to the reference genome, quantify transcript abundance, identify RNA Mod. and LincRNA
# ############################################################################################################################################################################################################################

# Paired end reads
if [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ]; then
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
      if [[ "$extension" =~ "fq.gz" ]]; then
        paired_fq_gz
      elif [[ "$extension" =~ "fastq.gz" ]]; then
        paired_fastq_gz
      elif [[ "$extension" =~ "fq" ]]; then
        echo "gzip" "$f"
        paired_fq
      elif [[ "$extension" =~ "fastq" ]]; then
        echo "gzip" "$f"
        paired_fastq
      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done
 fi

######### End ########