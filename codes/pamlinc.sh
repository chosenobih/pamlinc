#!/bin/bash
#Chosen Obih
#Script to process and analyze RNA-Seq data for transcript abundance, modified RNA and lincRNA

usage() {
      echo ""
      echo "Usage : sh $0 -g <reference_genome>  -a <reference_annotation> -l lib_type {-1 <left_reads> -2 <right_reads> | -u <single_reads> | -S <sra_id>} -o <output_folder for pipeline files> -p num_threads -t tophat -s star -q transcript_abundance_quantification -e evolinc_i -m HAMR"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########
  -g <reference genome fasta file>
  -a <reference genome annotation>
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
  -s star mapping #deactivates HAMR
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
transcript_abun_quant=0

while getopts ":g:a:l:1:2:u:o:S:p:htsqem:y:" opt; do
  case $opt in
    g)
    referencegenome=$OPTARG # Reference genome file
     ;;
    a)
    referenceannotation=$OPTARG # Reference genome annotation
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


paired_fastq_gz()
{
 filename=$(basename "$f" ".fq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
      mkdir $pipeline_output
        DIRECTORY=$pipeline_output/trimmomatic_out
          if [ ! -d "$DIRECTORY" ]; then
            mkdir $DIRECTORY
        echo "####################"
        echo "Running trimommatic"
        echo "####################"
  if [ "$seq_type" == "PE" ]; then
    echo "trimmomatic PE -thread $num_threads ${filename} ${filename2} ${filename}_forward_unpaired.fq.gz ${filename}_forward_paired.fq.gz ${filename2}_reverse_unpaired.fq.gz ${filename2}_reverse_paired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3
TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    trimmomatic PE -thread $num_threads ${filename} ${filename2} ${filename}_forward_unpaired.fq.gz ${filename}_forward_paired.fq.gz ${filename2}_reverse_unpaired.fq.gz ${filename2}_reverse_paired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3
TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    mv *_forward_paired *_forward_unpaired *_reverse_paired *_reverse_unpaired $DIRECTORY
        fi
  fi
}

# ############################################################################################################################################################################################################################
# # Reference genome/Index ###
# ############################################################################################################################################################################################################################

if [ "$tophat" != 0 ] && [ "$star" == 0 ]; then
  if [ ! -z "$index_folder" ]; then
    for i in $index_folder/*; do
        cp $i .
        fbname=$(basename "$i" .bt2 | cut -d. -f1)
    done
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
    echo "##########################################"
    echo "Building reference genome index for Tophat"
    echo "##########################################"
    echo "bowtie2-build -f $referencegenome ref_genome --threads $num_threads"
    bowtie2-build -f $referencegenome ref_genome --threads $num_threads
    echo "fbname=$(basename "ref_genome" .bt2 | cut -d. -f1)"
    fbname=$(basename "ref_genome" .bt2 | cut -d. -f1)
  fi
elif [ "$tophat" == 0 ] && [ "$star" != 0 ]; then
  if [ ! -z "$index_folder" ]; then
    for i in $index_folder/*; do
        cp $i -r .
    done
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
    echo "########################################"
    echo "Building reference genome index for STAR"
    echo "########################################"
    echo "STAR --runThreadN $num_threads --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles $referencegenome --sjdbGTFfile $referenceannotation --sjdbOverhang 100 --genomeSAindexNbases 13"
    STAR --runThreadN $num_threads --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles $referencegenome --sjdbGTFfile $referenceannotation --sjdbOverhang 100 --genomeSAindexNbases 13
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
   




######### End ########