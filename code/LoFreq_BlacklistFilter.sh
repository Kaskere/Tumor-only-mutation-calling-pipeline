#!/bin/bash
module load java
snpEff="/net/nfs/PAT/home/matias/tools/snpEff_4.3/"
#snpEff="/net/nfs/PAT/home/matias/tools/snpEff_4.1/"
SnpSift="/net/nfs/PAT/home/matias/tools/snpEff_4.3/SnpSift.jar"
#SnpSift="/net/nfs/PAT/home/matias/tools/snpEff_4.1/SnpSift.jar"
vcflib="/net/nfs/PAT/home/stef/tools/vcflib/bin/"
ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
manifest="/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/"
BED_blacklist="/net/nfs/PAT/home/matias/data/blacklist/PON_LoFreq_N24_min2POS_bq20.bed"
Gene_blacklist="/net/nfs/PAT/home/matias/data/blacklist/hypervariable_gene_blacklist.bed"

bgzip="/net/nfs/PAT/home/matias/tools/tabix-0.2.6/bgzip"
tabix="/net/nfs/PAT/home/matias/tools/tabix-0.2.6/tabix"


# Separate SNPs and Somatic variants from VCF file:

# create output dirs
if [ ! -d "vcf/blacklisted" ]; then mkdir "vcf/blacklisted" ; fi
if [ ! -d "vcf/not_blacklisted" ]; then mkdir "vcf/not_blacklisted" ; fi

for i in vcf/raw/*_lofreq.vcf
do 
	
	# in- and output files:
	bn=`basename $i`
	sname=${bn/_lofreq.vcf/}

	blacklisted_vcf="vcf/blacklisted/"$sname"_blacklisted.vcf"
	blacklisted_vcf_gz="vcf/blacklisted/"$sname"_blacklisted.vcf.gz"
	blacklisted_csv="vcf/blacklisted/"$sname"_blacklisted.csv"

	not_blacklisted_vcf="vcf/not_blacklisted/"$sname"_not_blacklisted.vcf"
	not_blacklisted_vcf_gz="vcf/not_blacklisted/"$sname"_not_blacklisted.vcf.gz"
	not_blacklisted_csv="vcf/not_blacklisted/"$sname"_not_blacklisted.csv"
	
	### 1 ### Separate blacklisted (Panel of normal) from non-PON:
	java -jar $SnpSift intervals -i $i $BED_blacklist $Gene_blacklist > $blacklisted_vcf
	java -jar $SnpSift intervals -i $i -x $BED_blacklist $Gene_blacklist > $not_blacklisted_vcf

	### 2 ### Extract desired columns from the vcf files (see Annotation_v1 script):
	java -jar $SnpSift extractFields -e "." $blacklisted_vcf \
	CHROM POS REF ALT DP AF > $blacklisted_csv

	java -jar $SnpSift extractFields -e "." $not_blacklisted_vcf \
	CHROM POS REF ALT DP AF > $not_blacklisted_csv

	# remove tmp files
	# rm $vcf_tmp

	# index the remaining vcf files:
	$bgzip -c $blacklisted_vcf > $blacklisted_vcf_gz
	$tabix -s1 -b2 -e2 $blacklisted_vcf_gz
	$bgzip -c $not_blacklisted_vcf > $not_blacklisted_vcf_gz
	$tabix -s1 -b2 -e2 $not_blacklisted_vcf_gz


done
# EOF


# Testing:
#i='raw/TAM15-34226FF_S21_L004_lofreq.vcf'

#java -jar $SnpSift intervals -i $i -x $BED_blacklist $Gene_blacklist > TAM15-34226FF_not_blacklisted.vcf 


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


















