## Workflow

1. FastQC – raw read quality assessment
2. Downsampling using seqtk
3. Adapter trimming using Trim Galore
4. Bisulfite-aware alignment using bwa-meth
5. BAM filtering and sorting using Sambamba
6. Duplicate marking using Picard
7. QC using CollectHsMetrics (how well the reads cover the Twist Panel target regions) 
8. Additional alignment QC using CollectMultipleMetrics
9. M-bias assessment using MethylDackel
10. Methylation extraction using MethylDackel 
