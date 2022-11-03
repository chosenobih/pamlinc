#!/bin/bash
#Chosen Obih
#I am using this script to keep record of the different tested sections of pamlinc as I work on them.






###################################################################################################################
# # Move files into output directory
###################################################################################################################

#### Trimmomatic
        DIRECTORY=$pipeline_output/trimmomatic_out
        if [ ! -d "$DIRECTORY" ]; then
            mkdir $DIRECTORY
          mv *_1P* *_1U* *_2P* *_2U* *trimlog* $DIRECTORY
	fi


############################################################################################################################################################################################################################
# # Trimmming, Mapping and greppiung unique reads ############################################################################################################################################################################################################################

# Paired end reads

paired_fastq_gz()
{
    filename=$(basename "$f" ".fastq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')

        echo "####################"
        echo "Running trimommatic"
        echo "####################"
          if [ "$seq_type" == "PE" ]; then
          echo "trimmomatic PE -threads $num_threads -trimlog ${filename3}_trimlog.txt ${filename} ${filename2} ${filename3}_1P.fastq.gz ${filename3}_1U.fastq.gz ${filename3}_2P.fastq.gz ${filename3}_2U.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
          trimmomatic PE -threads $num_threads -trimlog ${filename3}_trimlog.txt ${filename}.fastq.gz ${filename2}.fastq.gz ${filename3}_1P.fastq.gz ${filename3}_1U.fastq.gz ${filename3}_2P.fastq.gz ${filename3}_2U.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
          fi
}


 ############################################################################################################################################################################################################################
# # Check user supplied index folder for salmon, bowtie2 and star indexes. Build the indexes if not supplied ############################################################################################################################################################################################################################


if [ "$tophat" != 0 ] && [ "$star" == 0 ]; then
  if [ ! -z "$index_folder" ]; then
    for i in $index_folder/*; do
        mv $i -f .
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
        mv $i -f .
    done
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
    for i in $index_folder/*; do
        mv $i -f .
        fbname=$(basename "$i" . | cut -d. f1)
    done
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
  echo "####################################"
  echo "Building reference genome for salmon"
  echo "####################################"
  echo "salmon index -i salmon_index -t $referencegenome"
  salmon index -i salmon_index -t $referencegenome
fi

###################################################################################################
