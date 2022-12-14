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
# # Trimmming, Mapping and greppiung unique reads
############################################################################################################################################################################################################################

# Paired end reads

paired_fastq_gz()
{
    filename=$(basename "$f" ".fastq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')

          echo "######################"
          echo "Running trimommatic"
          echo "######################"
          
          if [ "$seq_type" == "PE" ]; then
          echo "trimmomatic PE -threads $num_threads -trimlog ${filename3}_trimlog.txt ${filename} ${filename2} ${filename3}_1P.fastq.gz ${filename3}_1U.fastq.gz ${filename3}_2P.fastq.gz ${filename3}_2U.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
          trimmomatic PE -threads $num_threads -trimlog ${filename3}_trimlog.txt ${filename}.fastq.gz ${filename2}.fastq.gz ${filename3}_1P.fastq.gz ${filename3}_1U.fastq.gz ${filename3}_2P.fastq.gz ${filename3}_2U.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
          fi
          
          if [ "$tophat" != 0 ] && [ "$star" == 0 ]; then
          echo "###################################"
          echo "Running tophat2 in paired end mode"
          echo "###################################"
          
          if [ "$seq_type" == "PE" ]; then
          echo "tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_fwd_tophat -G $referenceannotation ref_genome ${filename3}_1P.fastq.gz,${filename3}_1U.fastq.gz"
          tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_fwd_tophat -G $referenceannotation sbicolor ${filename3}_1P.fastq.gz,${filename3}_1U.fastq.gz
          echo "tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_rev_tophat -G $referenceannotation ref_genome ${filename3}_2P.fastq.gz,${filename3}_2U.fastq.gz"
          tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_rev_tophat -G $referenceannotation sbicolor ${filename3}_2P.fastq.gz,${filename3}_2U.fastq.gz
         
        
          echo "########################"
          echo "Converting .bam to .sam"
          echo "########################"
          echo "samtools view -h -@ $num_threads -o ${filename3}_fwd.sam ${filename3}_fwd/accepted_hits.bam"
          samtools view -h -@ $num_threads -o ${filename3}_fwd.sam ${filename3}_fwd/accepted_hits.bam
          echo "samtools view -h -o ${filename3}_rev.sam ${filename3}_rev/accepted_hits.bam"
          samtools view -h -o ${filename3}_rev.sam ${filename3}_rev/accepted_hits.bam
         
          echo "#######################"
          echo "Grepping unique reads"
          echo "#######################"
          echo "grep -P '^\@|NH:i:1$' ${filename3}_fwd.sam > ${filename3}_fwd_unique.sam"
          grep -P '^\@|NH:i:1$' ${filename3}_fwd.sam > ${filename3}_fwd_unique.sam
          echo "grep -P '^\@|NH:i:1$' ${filename3}_rev.sam > ${filename3}_rev_unique.sam"
          grep -P '^\@|NH:i:1$' ${filename3}_rev.sam > ${filename3}_rev_unique.sam
         
          echo "######################################################"
          echo "Converting .sam to .bam before running samtools sort"
          echo "######################################################"
          echo "samtools view -bSh -@ $num_threads ${filename3}_fwd_unique.sam > ${filename3}_fwd_unique.bam"
          samtools view -bSh -@ $num_threads ${filename3}_fwd_unique.sam > ${filename3}_fwd_unique.bam
          echo "samtools view -bSh -@ $num_threads ${filename3}_rev_unique.sam > ${filename3}_rev_unique.bam"
          samtools view -bSh -@ $num_threads ${filename3}_rev_unique.sam > ${filename3}_rev_unique.bam
         
          echo "#######################"
          echo "Sorting unique reads"
          echo "#######################"
          echo "samtools sort -@ $num_threads ${filename3}_fwd_unique.bam > ${filename3}_fwd_sorted.bam"
          samtools sort -@ $num_threads ${filename3}_fwd_unique.bam > ${filename3}_fwd_sorted.bam
          echo "samtools sort -@ $num_threads ${filename3}_rev_unique.bam > ${filename3}_rev_sorted.bam"
          samtools sort -@ $num_threads ${filename3}_rev_unique.bam > ${filename3}_rev_sorted.bam
         
          echo "########################"
          echo "Merging fwd and rev reads"
          echo "########################"
          echo "samtools merge -@ $num_threads -f ${filename3}_merged.bam ${filename3}_fwd_sorted.bam ${filename3}_rev_sorted.bam"
          samtools merge -@ $num_threads -f ${filename3}_merged.bam ${filename3}_fwd_sorted.bam ${filename3}_rev_sorted.bam

	  echo "######################################################"
          echo "Resolving spliced alignments"
          echo "######################################################"
          echo "picard AddOrReplaceReadGroups I=${filename3}_merged.bam O=${filename3}_RG.bam ID=${filename3} LB=D4 PL=illumina PU=HWUSI-EAS1814:28:2 SM=${filename3}"
          picard AddOrReplaceReadGroups I=${filename3}_merged.bam O=${filename3}_RG.bam ID=${filename3} LB=D4 PL=illumina PU=HWUSI-EAS1814:28:2 SM=${filename3}
          echo "picard ReorderSam I=${filename3}_RG.bam O=${filename3}_RGO.bam R=$referencegenome"
          picard ReorderSam I=${filename3}_RG.bam O=${filename3}_RGO.bam R=$referencegenome
          echo "samtools index ${filename3}_RGO.bam ${filename3}_RGO.bam.bai"
          samtools index ${filename3}_RGO.bam ${filename3}_RGO.bam.bai
          echo "java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R $referencegenome -I ${filename3}_RGO.bam -o ${filename3}_resolvedalig.bam -U ALLOW_N_CIGAR_READS"
          java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R $referencegenome -I ${filename3}_RGO.bam -o ${filename3}_resolvedalig.bam -U ALLOW_N_CIGAR_READS
	  
	  echo "######################################################"
          echo "Running HAMR"
          echo "######################################################"
	  echo "singularity run --cleanenv hamr_xi_1.4.sif -fe ${filename3}_resolvedalig.bam $referencegenome hamr_model/euk_trna_mods.Rdata ${filename3}_HAMR ${filename3} 30 10 0.01 H4 1 .05 .05"
	  singularity run --cleanenv hamr_xi_1.4.sif -fe ${filename3}_resolvedalig.bam $referencegenome hamr_model/euk_trna_mods.Rdata ${filename3}_HAMR ${filename3} 30 10 0.01 H4 1 .05 .05
		
}


 ############################################################################################################################################################################################################################
# # Check if user supplied index folder for salmon, bowtie2 and star indexes. Build the indexes if not supplied.
#############################################################################################################################################################################################################################

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
 ############################################################################################################################################################################################################################
# # Check that the input fastq files has the appropriate extension and then trim reads, align the reads to the reference genome, quantify transcript abundance, identify RNA Mod. and LincRNA
#############################################################################################################################################################################################################################

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
