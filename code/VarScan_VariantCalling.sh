#!/bin/bash
module load java
module load samtools

VarScan="/net/nfs/PAT/home/matias/tools/VarScan/VarScan.v2.3.7.jar"
#VarScan="/net/nfs/PAT/lib/VarScan/VarScan.v2.2.11.jar"
vcflib="/net/nfs/PAT/home/stef/tools/vcflib/bin/"
ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
#manifest="/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/"


# create vcf dir
if [ ! -d "vcf" ]; then mkdir "vcf" ; fi	
if [ ! -d "vcf/raw" ]; then mkdir "vcf/raw" ; fi	

#for i in bam/*_coordsorted.bam
#for i in bam/*_recal.bam
for i in bam/lofreq_indelq/*_indelq_sorted.bam

do
        # in- and output files:
        bn=`basename $i`
		sname=${bn/_indelq_sorted.bam/}

        raw_snps="vcf/"$sname"_varscan_snps.vcf"
        raw_indels="vcf/"$sname"_varscan_indels.vcf"
        tmp="vcf/"$sname"_varscan_tmp.vcf"
        out="vcf/raw/"$sname"_varscan.vcf"

        # call SNPs
        echo `date` " [$$] Calling variants (VarScan2) for: " $i
        samtools mpileup -f $ref $i | java -jar $VarScan mpileup2snp \
                --min-coverage 20 \
                --min-reads2 6 \
                --min-avg-qual 20 \
                --min-var-freq 0.05 \
                --p-value 0.01 \
                --strand-filter 1 \
                --output-vcf 1 \
                --variants 0 \
                > $raw_snps

        # call INDELs
        samtools mpileup -f $ref $i | java -jar $VarScan mpileup2indel \
                --min-coverage 20 \
                --min-reads2 6 \
                --min-avg-qual 20 \
                --min-var-freq 0.05 \
                --p-value 0.01 \
                --strand-filter 1 \
                --output-vcf 1 \
                --variants 0 \
		> $raw_indels
		
		# merge raw snp and indel vcfs
		$vcflib/vcfcombine $raw_snps $raw_indels > $tmp
		# Remove '%' from FREQ format-field and rename Integer > Float for FREQ-field:
		sed 's/\%//' $tmp | sed 's/FREQ\,Number\=1\,Type\=String/FREQ\,Number\=1\,Type\=Float/' > $out
		# remove redundant files:
		rm $tmp $raw_snps $raw_indels

done


# --min-avg-qual  Minimum base quality at a position to count a read [15]

# compared to Lofreq, settings that are similar:
	# min-cov (20)
	# min-base qual (20)
	# min variant supporting reads: 6 (>3 alt-forward & >3 alt-rev)
	# min-vaf (0.05)
	# sign p-value (0.01) 
	# strand filter (VarScan 1; LoFreq sb-alpha 0.05)

# compared to Lofreq, settings that are not similar/not used:
	# min-mapping qual (30)
	# max-depth (1000)






