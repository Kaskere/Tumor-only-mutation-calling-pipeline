#!/bin/bash

##############################################################################################
# Select aligned reads from bam file that overlap with capture target regions 
# for quicker processing (for somatic mutation analysis only)
# downside: SNPs outside of target regions (introns) are not called/evaluated
# if not used, still after each calling tool, only mutations in regions of interest are kept.

##############################################################################################
# DEPENDENCIES:
# 

module load java
module load samtools

picard="/net/nfs/PAT/lib/picard-tools/picard-tools-1.61/"

# mutation target regions:
target_bed='/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/BCNHLv2_allExons.bed'

##############################################################################################
# TODO: bam directory/ *{string_bamfile_extension}
# TODO: use the same *{string_bamfile_extension}
# TODO: use bam directory from input

for i in bam/*_recal.bam 					
do
# in- and output files:
bn=`basename $i`
sname=${bn/_recal.bam/}
bamout="bam/"$sname"_target_exons.bam"

# Select aligned reads from bam file that overlap with target regions 
samtools view -L $target_bed $i -b -o $bamout
# create bam index
java -jar $picard/BuildBamIndex.jar I=$bamout

done
# EOF




##############################################################################################
# TEST
##############################################################################################

# mutation target regions:
#target_bed='/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/BCNHLv2_allExons.bed'
#bamfile='../bam/TAM15-34226_S22_L004_recal.bam'
#samtools view -L $target_bed $bamfile -b -o real.bam
# check file length and size
#samtools view $bamfile | wc -l # 56.154.131, size: 6.5 Gb
#samtools view real.bam | wc -l # 9.235.023 , size: 744M
#less out.bam | wc -l # (uncompressed, when option -b is not used) 9.235.023. size: 3.4G


