#!/bin/bash
#Chosen Obih
#Script to process and analyze RNA-Seq data for transcript abundance, modified RNA and lincRNA

usage() {
      echo ""
      echo "Usage : sh $0 -g <reference_genome>  -a <reference_annotation> -i <index_folder> -l lib_type {-1 <left_reads> -2 <right_reads> | -u <single_reads> | -S <sra_id>} -o <output_folder for pipeline files> -p num_threads -t tophat -s star -q transcript_abundance_quantification -e evolinc_i -m HAMR"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########
  -g <reference genome fasta file>
  -a <reference genome annotation>
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
  -y type of reads (single end or paired end) #denoted as "SE" or "PE", include double quotation on command line
  -m HAMR
  -e evolinc_i
EOF
    exit 0
}

star=0
tophat=0
referencegenome=0
referenceannotation=0

while getopts ":g:a:i:l:1:2:u:o:S:p:htsqem:y:" opt; do
  case $opt in
    g)
    referencegenome=$OPTARG # Reference genome file
     ;;
    a)
    referenceannotation=$OPTARG # Reference genome annotation
     ;;
    i)
    index_folder=$OPTARG # Input folder
     ;;
    l)
    lib_type=$OPTARG # Library type
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
    t)
    tophat=$OPTARG # tophat
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
mkdir $pipeline_output

paired_fastq_gz()
{
 filename=$(basename "$f" ".fastq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
        DIRECTORY=$pipeline_output/trimmomatic_out
          if [ ! -d "$DIRECTORY" ]; then
            mkdir $DIRECTORY
        echo "####################"
        echo "Running trimommatic"
        echo "####################"
  if [ "$seq_type" == "PE" ]; then
    echo "trimmomatic PE -threads $num_threads -trimlog ${filename3}_trimlog.txt ${filename} ${filename2} ${filename3}_1P.fq.gz ${filename3}_1U.fq.gz ${filename3}_2P.fq.gz ${filename3}_2U.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    trimmomatic PE -threads $num_threads -trimlog ${filename3}_trimlog.txt ${filename} ${filename2} ${filename3}_1P.fq.gz ${filename3}_1U.fq.gz ${filename3}_2P.fq.gz ${filename3}_2U.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    mv *_1P* *_1U* *_2P* *_2U* *trimlog* $DIRECTORY
        fi
  fi
}


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