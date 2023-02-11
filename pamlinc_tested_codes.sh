#!/bin/bash
#Chosen Obih
#Script to process and analyze RNA-Seq data for transcript abundance, modified RNA and lincRNA

usage() {
      echo ""
      echo "Usage : sh $0 -g <reference_genome>  -a <reference_annotation> -A <reference_annotation type> -i <index_folder> -l lib_type {-1 <left_reads> -2 <right_reads> | -u <single_reads> | -S <sra_id>} -o <output_folder for pipeline files> -p num_threads -d reads_mismatches -t tophat -s star -q transcript_abundance_quantification -e evolinc_i -m HAMR -r <gene_attribute> -n <strandedness> -k <feature_type>"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########
  -g <reference genome fasta file>
  -a <reference genome annotation>
  -A <ref_annot type> ("GTF" or "GFF3" supported include double quotation on command line)
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
  -b reads_mismatches (% reads mismatches to allow. Needed for tophat2)
  -m HAMR
  -e evolinc_i
  -k feature_type #Feature type (Default is exon)
  -r gene attribute (Default is gene_id)
  -n strandedness (Default is 0 (unstranded), 1 (stranded), 2 (reversely stranded)
  -z type of container (docker or singularity) #denoted as "D" or "S", include double quotation on command line
EOF
    exit 0
}

star=0
tophat=0
referencegenome=0
referenceannotation=0

while getopts ":g:a:i:l:1:2:u:o:S:p:d:htsqem:k:r:n:y:" opt; do
  case $opt in
    g)
    referencegenome=$OPTARG # Reference genome file
     ;;
    a)
    referenceannotation=$OPTARG # Reference genome annotation
     ;;
    A)
    ref_annot_type=$OPTARG # Reference genome annotation type (GTF or GFF3 supported)
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
     z)
    container_type=$OPTARG # Type of container ("D" or "S". Option is required when running pamlinc as a container)
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


###################################################################################################################
# # Check reference genome annotation file type and convert to .gtf if .gff file was supplied by user.
###################################################################################################################
if (grep -q -E 'transcript_id | gene_id' $referenceannotation); then
    echo "$referenceannotation is in .gtf format"
    else
    gffread $referenceannotation -T -o ref_annotation.gtf
    referenceannotation=ref_annotation.gtf
fi

###################################################################################################################
# # Pull required images for HAMR, EVOLINC-I and GATK
###################################################################################################################
if [ "$container_type" == "D" ]; then
echo "###################################################"
echo "Pulling docker images for HAMR, EVOLINC-I and GATK"
echo "###################################################"
    docker pull reetututeja/hamr_xi:1.4
    docker pull evolinc/evolinc-i:1.7.5
    docker pull broadinstitute/gatk3:3.5-0
elif [ "$container_type" == "S" ]; then
echo "########################################################"
echo "Pulling singularity images for HAMR, EVOLINC-I and GATK"
echo "########################################################"
    singularity pull docker://reetututeja/hamr_xi:1.4
    singularity pull docker://evolinc/evolinc-i:1.7.5
    singularity pull docker://broadinstitute/gatk3:3.5-0
fi

###################################################################################################################
# # House keeping - move files into output directory
###################################################################################################################

#### Trimmomatic
#        DIRECTORY=$pipeline_output/trimmomatic_out
#        if [ ! -d "$DIRECTORY" ]; then
#            mkdir $DIRECTORY
#          mv *_1P* *_1U* *_2P* *_2U* *trimlog* $DIRECTORY
#	fi
	
##################################################################################################################
# # Transcript abundance quantification
###################################################################################################################

transcript_quantification()
{
      if [ "$transcript_abun_quant" != 0 ]; then
      echo "###########################################################################"
      echo "Running salmon and featureCounts to quantify transcript"
      echo "###########################################################################"    
        if [ "$seq_type" == "PE" ]; then
        echo "salmon quant -l A -a ${filename3}_merged.bam -o ${filename3}_salmon -t $referencegenome -p $num_threads"
        salmon quant -l A -a ${filename3}_merged.bam -o ${filename3}_salmon -t $referencegenome -p $num_threads
        echo "featureCounts -p -T $num_threads -t $feature_type -g $gene_attribute -s $strandedness -a $referenceannotation -o ${filename3}_featurecount.txt ${filename3}_merged.bam"
        featureCounts -p -T $num_threads -t $feature_type -g $gene_attribute -s $strandedness -a $referenceannotation -o ${filename3}_featurecount.txt ${filename3}_merged.bam
        #elif [ "$seq_type" == "SE" ]; then
        #echo "salmon quant -l A -a ${filename3}.bam -o ${filename3}_salmon -t $referencegenome -p $num_threads"
        #salmon quant -l A -a ${filename3}.bam -o ${filename3}_salmon -t $referencegenome -p $num_threads
        #echo "featureCounts -T $num_threads -t $feature_type -g $gene_attribute -s $strandedness -a $referenceannotation -o ${filename3}_featurecount.txt ${filename3}_merged.bam"
        #featureCounts -T $num_threads -t $feature_type -g $gene_attribute -s $strandedness -a $referenceannotation -o ${filename3}_featurecount.txt ${filename3}_merged.bam
        fi 
      fi
}

##################################################################################################################
# # lincRNA identification
###################################################################################################################

 lincRNA_annotation()
{
      if [ "$evolinc_i" != 0 ]; then
      echo "###########################################################################"
      echo "Converting .bam file(s) containing uniquely mapped and sorted reads to .gtf"
      echo "###########################################################################"
        
        if [ "$lib_type" == fr-secondstrand ]; then      
        echo "stringtie ${filename3}_merged.bam -o ${filename3}_merged.gtf -G $referenceannotation -p $num_threads --fr"
        stringtie ${filename3}_merged.bam -o ${filename3}_merged.gtf -G $referenceannotation -p $num_threads --fr
        echo "cuffcompare ${filename3}_merged.gtf -r $referenceannotation -s $referencegenome -T -o ${filename3}"
        cuffcompare ${filename3}_merged.gtf -r $referenceannotation -s $referencegenome -T -o ${filename3}
        echo "singularity run -B $(pwd):/mnt --pwd /mnt evolinc-i_1.7.5.sif -c ./${filename3}.combined.gtf -g ./$referencegenome -u ./$referenceannotation -r ./$referenceannotation -n $num_threads -o ./${filename3}_lincRNA"
        singularity run -B $(pwd):/mnt --pwd /mnt evolinc-i_1.7.5.sif -c ./${filename3}.combined.gtf -g ./$referencegenome -u ./$referenceannotation -r ./$referenceannotation -n $num_threads -o ./${filename3}_lincRNA
        elif [ "$lib_type" == fr-firststrand ]; then
        echo "stringtie ${filename3}_merged.bam -o ${filename3}_merged.gtf -G $referenceannotation -p $num_threads --rf"
        stringtie ${filename3}_merged.bam -o ${filename3}_merged.gtf -G $referenceannotation -p $num_threads --rf
        echo "cuffcompare ${filename3}_merged.gtf -r $referenceannotation -s $referencegenome -T -o ${filename3}"
        cuffcompare ${filename3}_merged.gtf -r $referenceannotation -s $referencegenome -T -o ${filename3}
        echo "singularity run -B $(pwd):/mnt --pwd /mnt evolinc-i_1.7.5.sif -c ./${filename3}.combined.gtf -g ./$referencegenome -u ./$referenceannotation -r ./$referenceannotation -n $num_threads -o ./${filename3}_lincRNA"
        singularity run -B $(pwd):/mnt --pwd /mnt evolinc-i_1.7.5.sif -c ./${filename3}.combined.gtf -g ./$referencegenome -u ./$referenceannotation -r ./$referenceannotation -n $num_threads -o ./${filename3}_lincRNA
        fi
      fi
}

############################################################################################################################################################################################################################
# # Trimmming, Mapping and greppiung unique reads 
############################################################################################################################################################################################################################

# Paired end reads

paired_fastq_gz()
{
    filename=$(basename "$f" ".fastq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
          
          if [ "$seq_type" == "PE" ]; then
          echo "######################"
          echo "Trimming input read(s)"
          echo "######################"
          
          echo "trimmomatic PE -threads $num_threads -trimlog ${filename3}_trimlog.txt ${filename} ${filename2} ${filename3}_1P.fastq.gz ${filename3}_1U.fastq.gz ${filename3}_2P.fastq.gz ${filename3}_2U.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
          trimmomatic PE -threads $num_threads -trimlog ${filename3}_trimlog.txt ${filename}.fastq.gz ${filename2}.fastq.gz ${filename3}_1P.fastq.gz ${filename3}_1U.fastq.gz ${filename3}_2P.fastq.gz ${filename3}_2U.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
          fi
          
          if [ "$tophat" != 0 ] && [ "$star" == 0 ]; then
          echo "###################################"
          echo "Running tophat2 in paired end mode"
          echo "###################################"
          
          if [ "$seq_type" == "PE" ]; then
          echo "tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_fwd_tophat -G $referenceannotation sbicolor ${filename3}_1P.fastq.gz,${filename3}_1U.fastq.gz"
          tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_fwd_tophat -G $referenceannotation sbicolor ${filename3}_1P.fastq.gz,${filename3}_1U.fastq.gz
          echo "tophat2 -p $num_threads --library-type $lib_type --read-mismatches $reads_mismatches --read-edit-dist $reads_mismatches --max-multihits 10 --b2-very-sensitive --transcriptome-max-hits 10 --no-coverage-search --output-dir ${filename3}_rev_tophat -G $referenceannotation sbicolor ${filename3}_2P.fastq.gz,${filename3}_2U.fastq.gz"
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
          singularity run --cleanenv gatk3_3.5-0.sif java -Xmx8g -jar ./GenomeAnalysisTK.jar -T SplitNCigarReads -R $referencegenome -I ${filename3}_RGO.bam -o ${filename3}_resolvedalig.bam -U ALLOW_N_CIGAR_READS
	  
	        echo "######################################################"
          echo "Running HAMR"
          echo "######################################################"
	        echo "singularity run --cleanenv hamr_xi_1.4.sif -fe ${filename3}_resolvedalig.bam $referencegenome hamr_model/euk_trna_mods.Rdata ${filename3}_HAMR ${filename3} 30 10 0.01 H4 1 .05 .05"
	        singularity run --cleanenv hamr_xi_1.4.sif -fe ${filename3}_resolvedalig.bam $referencegenome hamr_model/euk_trna_mods.Rdata ${filename3}_HAMR ${filename3} 30 10 0.01 H4 1 .05 .05
	        lincRNA_annotation
	        transcript_quantification
          fi
          fi
}

#############################################################################################################################################################################################################################
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
    echo "bowtie2-build --threads -f $referencegenome ref_genome"
    bowtie2-build -f $referencegenome ref_genome
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
  fi
fi

if [ "$transcript_abun_quant" != 0 ]; then
  if [ ! -z "$index_folder" ]; then
    for i in $index_folder/*; do
        mv $i -f .
    done
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
  echo "####################################"
  echo "Building reference genome for salmon"
  echo "####################################"
  echo "salmon index -i salmon_index -t $referenceannotation"
  salmon index -i salmon_index -t $referenceannotation
  fi
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
######### End ########