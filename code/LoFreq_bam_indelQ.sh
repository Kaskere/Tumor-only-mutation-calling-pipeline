#!/bin/bash
module load java
module load samtools
ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
picard="/net/nfs/PAT/lib/picard-tools/picard-tools-1.61/"

# create vcf dir
if [ ! -d "bam/lofreq_indelq" ]; then mkdir "bam/lofreq_indelq" ; fi	

#for i in bam/*_recal.bam
for i in bam/*_target_exons.bam

do
# in- and output files:
bn=`basename $i`
#sname=${bn/_recal.bam/}
sname=${bn/_target_exons.bam/}
indelq_bam="bam/lofreq_indelq/"$sname"_indelq.bam"	
indelq_sorted_bam="bam/lofreq_indelq/"$sname"_indelq_sorted.bam"	

# lofreq indelqual --dindel 
echo `date` " [$$] Calculating indel qualities (LoFreqStar) for: " $i		
lofreq indelqual --dindel -f $ref $i -o $indelq_bam

# picard sort bam + create index				
java -Xmx4g -jar $picard/SortSam.jar I=$indelq_bam O=$indelq_sorted_bam SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR=`pwd`/tmp

# remove indelq bam
rm $indelq_bam
# rm $i #(target_exons.bam)

done


