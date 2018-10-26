#!/bin/bash
module load java
#module load samtools
#module load bwa

picard="/net/nfs/PAT/lib/picard-tools/picard-tools-1.61/"
lofreqPar='/net/nfs/PAT/home/matias/tools/miniconda2/bin/lofreq2_call_pparallel.py'
SnpSift="/net/nfs/PAT/home/matias/tools/snpEff_4.3/SnpSift.jar"
THREADS=8
ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
#ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
manifest="/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2"

dbSnp="/net/nfs/PAT/home/matias/dbSNP/b150/All_20170710.vcf.gz"

# create vcf dir
if [ ! -d "vcf" ]; then mkdir "vcf" ; fi	
if [ ! -d "vcf/raw" ]; then mkdir "vcf/raw" ; fi	

for i in bam/lofreq_indelq/*_indelq_sorted.bam

do
# in- and output files:
bn=`basename $i`
sname=${bn/_indelq_sorted.bam/}

raw_SNVs="vcf/raw/"$sname"_lofreq_raw.vcf"
tmp_vcf="vcf/raw/"$sname"_lofreq_tmp.vcf"
lofreq_vcf="vcf/raw/"$sname"_lofreq.vcf"
#echo $indelq_sorted_bam
#echo $raw_SNVs
#done

# lofreq call --call-indels
echo `date` " [$$] Calling variants (LoFreqStar) for: " $i		
python $lofreqPar --pp-threads $THREADS --call-indels --verbose $i -f $ref -o $raw_SNVs \
--min-cov 20 \
--min-mq 30 \
--min-bq 20 \
--min-alt-bq 20 \
--max-depth 1000 \
--sig 0.01

# -l $manifest/BCNHLv2_primary_coord.bed \ # the manifest file is not being used, since the bam files are already filtered for reads on target regions

# lofreq filter 
echo `date` " [$$] Filtering variants (LoFreqStar) for: " $i		
lofreq filter --verbose --cov-min 10 --af-min 0.05 --sb-alpha 0.05 --sb-incl-indels -i $raw_SNVs -o $tmp_vcf
# TODO: remove --cov-min 10 ; since it is already used in lofreq call

# Alternative bases / Strand Bias Filter
# Minimal 3 alt-forward and 3 alt-reverse bases">
java -jar $SnpSift filter -f $tmp_vcf "(DP4[2]>3) & (DP4[3]>3)" > $lofreq_vcf

# remove tmp files
rm $tmp_vcf

# TODO: add homopolymer filter (if necessary)
# java -jar $SnpSift filter -f $tmp_vcf "(HRUN> ...)" > $lofreq_vcf

done

# info from verbose filter
# --snvqual-thresh 77 --indelqual-thresh 61
# Setting default SB filtering method to FDR
# At least one type of multiple testing correction requested. Doing first pass of vcf
# MTC application completed
# Successful exit.




# --sb-alpha FLOAT          Multiple testing correction pvalue threshold


# FROM VCF
##FILTER=<ID=min_snvqual_63,Description="Minimum SNV Quality (Phred) 63">
##FILTER=<ID=min_indelqual_47,Description="Minimum Indel Quality (Phred) 47">

##FILTER=<ID=min_snvqual_77,Description="Minimum SNV Quality (Phred) 77">
##FILTER=<ID=min_indelqual_61,Description="Minimum Indel Quality (Phred) 61">

##FILTER=<ID=sb_fdr,Description="Strand-Bias Multiple Testing Correction: fdr corr. pvalue > 0.001000">






#i='/net/nfs/PAT/analysis/MPS-310/SNVanalysis_normals/bam/PMBC-1_S1_recal.bam'
# testing
#python $lofreqPar --pp-threads $THREADS --call-indels --verbose $i -f $ref -o $raw_SNVs \
#-l $manifest/BCNHLv2_primary_coord.bed \
#--min-cov 10 \
#--max-depth 1000 \
#--sig 0.05 

#lofreq call --call-indels --verbose $i -f $ref -o $raw_SNVs \
#-l $manifest/BCNHLv2_primary_coord.bed \
#--min-cov 10 \
#--max-depth 1000 \
#--sig 0.05 

