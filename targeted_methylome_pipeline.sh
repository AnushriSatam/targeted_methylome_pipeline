#!/bin/bash
set -e
set -o pipefail

if [ $# -ne 2 ];then
echo "usage: $0 <R1> <R2>"
exit 1
fi

R1=$1
R2=$2
BASENAME=$(basename "$R1")
READS=101790680 #required for the downsampling step. This is the standardised read count used for downsampling 

#genomic location
genome=hg19_methylom_index/hg19_exome.fasta
interval=hg19_methylom_index/Twist_Human_Methylome_Panel_hg19.intervals
#software location
#the location must be mentioned here.
seqtk=""
trim_galore="" 
bwa_meth=""
sambamba=""
picard="" 
MethylDackel=""
#QC
echo "qc"
fastqc --threads 16 $R1 $R2

#DOWNSAMPLE
echo "downsampling"
$seqtk sample -2 -s 42 "$R1" $READS | gzip > ${BASENAME}_R1_down.fastq.gz
$seqtk sample -2 -s 42 "$R2" $READS | gzip > ${BASENAME}_R2_down.fastq.gz

#trimming
#2Color handles two colour artifacts
echo "trimming"
$trim_galore \
  --paired \
  --gzip \
  --cores 16 \
  --2colour 20 \ 
  ${BASENAME}_R1_down.fastq.gz \
  ${BASENAME}_R2_down.fastq.gz

echo "aligning"
#bwa would not be able to handle the bisulfite conversion
$bwa_meth \   
  --reference $genome \
  --read-group "@RG\tID:${BASENAME}\tPL:illumina\tLB:${BASENAME}\tSM:${BASENAME}" \
  --threads 40 \
  ${BASENAME}_R1_down_val_1.fq.gz \
  ${BASENAME}_R2_down_val_2.fq.gz | \
$sambamba view -S -f bam -t 16 --filter "not secondary_alignment and not supplementary and proper_pair and mapping_quality > 0" /dev/stdin | \ 
$sambamba sort -t 16 -o ${BASENAME}_sorted.bam /dev/stdin

#markduplicates 
$picard MarkDuplicates \
  -I ${BASENAME}_sorted.bam \
  -O ${BASENAME}_sorted.markdup.bam \
  -M ${BASENAME}_markdup_metrics.txt \
  -R $genome \
  --ASSUME_SORT_ORDER coordinate \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --MAX_RECORDS_IN_RAM 1000 \
  --CREATE_INDEX false

samtools index -@ 16 ${BASENAME}_sorted.markdup.bam


#Collect HsMetrics-how well did the reads hit the intended target regions (twist panel) 
$picard CollectHsMetrics \
-I ${BASENAME}_sorted.markdup.bam \
-O ${BASENAME}_metrics.txt \
-R $genome \
--BAIT_INTERVALS $interval \
--TARGET_INTERVALS $interval \
--MINIMUM_MAPPING_QUALITY 20 \
--NEAR_DISTANCE 500 \
--COVERAGE_CAP 1000 \
--PER_TARGET_COVERAGE per_target.txt

#additional metrices
$picard CollectMultipleMetrics \
-I ${BASENAME}_sorted.markdup.bam \
-O ${BASENAME}_metrics \
-R $genome \
--PROGRAM null \
--PROGRAM CollectGcBiasMetrics \
--PROGRAM CollectInsertSizeMetrics \
--PROGRAM CollectAlignmentSummaryMetrics



BAM="${BASENAME}_sorted.markdup.bam"

#checking for methylation bias.(decreased or increased methylation is  observed at the ends of the reads)
$MethylDackel mbias $genome $BAM ${BASENAME}_methylome
#however the values used in the below step is based on Twist internal dataset. 

$MethylDackel extract \
--minDepth 10 \    
--maxVariantFrac 0.25 \
--OT 0,0,0,98 \
--OB 0,0,3,0 \
--mergeContext \
$genome \
$BAM -o ${BASENAME}_methylome


echo "pipeline completed."

