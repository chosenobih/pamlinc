#!/bin/bash
#Chosen Obih
#I am using this script to keep record of the different sections of pamlinc as I work on them.


# ############################################################################################################################################################################################################################
# # Check user supplied index folder for salmon, bowtie2 and star indexes. Build the indexes if not supplied
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

if [ "$transcript_abun_quant" != 0 ]; then
  if [ ! -z "$index_folder" ]; then
    for i in $index_folder/*; do
        cp $i -r .
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

