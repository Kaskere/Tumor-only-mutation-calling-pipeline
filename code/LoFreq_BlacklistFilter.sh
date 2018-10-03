#!/bin/bash
module load java
snpEff="/net/nfs/PAT/home/matias/tools/snpEff_4.3/"
#snpEff="/net/nfs/PAT/home/matias/tools/snpEff_4.1/"
SnpSift="/net/nfs/PAT/home/matias/tools/snpEff_4.3/SnpSift.jar"
#SnpSift="/net/nfs/PAT/home/matias/tools/snpEff_4.1/SnpSift.jar"
vcflib="/net/nfs/PAT/home/stef/tools/vcflib/bin/"
ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
manifest="/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/"

#PON_blacklist="/net/nfs/PAT/home/matias/data/blacklist/PON_LoFreq_N24.vcf.gz"
#BED_blacklist="/net/nfs/PAT/analysis/MPS-310/SNVanalysis_normals/LoFreq/PON_LoFreq/PON_LoFreq_N24.bed"
alt_BED_blacklist="/net/nfs/PAT/home/matias/data/blacklist/PON_LoFreq_N24_Adjusted_bq20.bed"
#alt_BED_blacklist="/net/nfs/PAT/analysis/MPS-310/SNVanalysis_normals/LoFreq/PON_LoFreq/PON_LoFreq_N24_ALT.bed"
# alt_BED_blacklist is the same as the BED_blacklist from Lofreq, except all positions are listed -1. 

# Separate SNPs and Somatic variants from VCF file:

#mkdir somatic
# create vcf dir
#if [ ! -d "vcf" ]; then mkdir "vcf" ; fi	
if [ ! -d "vcf/blacklisted" ]; then mkdir "vcf/blacklisted" ; fi
if [ ! -d "vcf/not_blacklisted" ]; then mkdir "vcf/not_blacklisted" ; fi

for i in vcf/raw/*_lofreq.vcf
do 
	
	# in- and output files:
	bn=`basename $i`
	sname=${bn/_lofreq.vcf/}

	#vcf_tmp="vcf/"$sname"_tmp.vcf"

	blacklisted_vcf="vcf/blacklisted/"$sname"_blacklisted.vcf"
	blacklisted_csv="vcf/blacklisted/"$sname"_blacklisted.csv"

	not_blacklisted_vcf="vcf/not_blacklisted/"$sname"_not_blacklisted.vcf"
	not_blacklisted_csv="vcf/not_blacklisted/"$sname"_not_blacklisted.csv"
	
	### 1 ### Separate blacklisted (Panel of normal) from non-PON:
	java -jar $SnpSift intervals -i $i $alt_BED_blacklist > $blacklisted_vcf
	java -jar $SnpSift intervals -i $i -x $alt_BED_blacklist > $not_blacklisted_vcf


	### 2 ### Extract desired columns from the vcf files (see Annotation_v1 script):
	java -jar $SnpSift extractFields -e "." $blacklisted_vcf \
	CHROM POS REF ALT DP AF > $blacklisted_csv

	java -jar $SnpSift extractFields -e "." $not_blacklisted_vcf \
	CHROM POS REF ALT DP AF > $not_blacklisted_csv

	# remove tmp files
	# rm $vcf_tmp
	#$vcf_tmp2 $vcf_tmp3

	# index the remaining vcf files:
	#$bgzip -c $vcf_out > $vcf_out_gz
	#$tabix -p vcf $vcf_out_gz
	
	
done
# EOF


# Testing:
#i='raw/TAM15-34226FF_S21_L004_lofreq.vcf'

# with the adjusted BED files (positions -1) ; it seems to give desirable results
#java -jar $SnpSift intervals -i $i $alt_BED_blacklist > TAM15-34226FF_blacklisted.vcf 
#java -jar $SnpSift intervals -i $i -x $alt_BED_blacklist > TAM15-34226FF_not_blacklisted.vcf 


# this seems to give the same output as : java -jar $SnpSift filter -f $i "(na PON_BLACKLIST)" > $not_blacklisted_vcf
#java -jar $SnpSift intervals -i $i $BED_blacklist > TAM15-34226FF_blacklisted.vcf 
#java -jar $SnpSift intervals -i $i -x $BED_blacklist > TAM15-34226FF_not_blacklisted.vcf 	



# try using bedtools
#bedtools='/net/nfs/PAT/home/matias/tools/bedtools2/bin'
#$bedtools/intersectBed -a input.vcf -b /path/to/my.interval.bed -header > output.vcf
#$bedtools/intersectBed -a $i -b $BED_blacklist -header > output.vcf
#$bedtools/intersectBed -a $i -b $alt_BED_blacklist -header > output.vcf


















