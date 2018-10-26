
#!/bin/bash
module load java
snpEff="/net/nfs/PAT/home/matias/tools/snpEff_4.3/"
SnpSift="/net/nfs/PAT/home/matias/tools/snpEff_4.3/SnpSift.jar"
vcflib="/net/nfs/PAT/home/stef/tools/vcflib/bin/"
ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
manifest="/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/"

# Separate SNPs and Somatic variants from VCF file:

# create vcf dir
if [ ! -d "vcf/snp" ]; then mkdir "vcf/snp" ; fi	
if [ ! -d "vcf/somatic" ]; then mkdir "vcf/somatic" ; fi	
if [ ! -d "vcf/somatic/all" ]; then mkdir "vcf/somatic/all" ; fi	
if [ ! -d "vcf/somatic/functional" ]; then mkdir "vcf/somatic/functional" ; fi	
if [ ! -d "vcf/somatic/non_functional" ]; then mkdir "vcf/somatic/non_functional" ; fi	

for i in vcf/annotated/*_annot.vcf
do 
	
	# in- and output files:
	bn=`basename $i`
	sname=${bn/_annot.vcf/}

	blacklisted_vcf="vcf/blacklisted/"$sname"_blacklisted.vcf"
	blacklisted_csv="vcf/blacklisted/"$sname"_blacklisted.csv"

	not_blacklisted_vcf="vcf/not_blacklisted/"$sname"_not_blacklisted.vcf"
	not_blacklisted_csv="vcf/not_blacklisted/"$sname"_not_blacklisted.csv"
	
	somatic_vcf="vcf/somatic/all/"$sname"_all_somatics.vcf"
	somatic_csv="vcf/somatic/all/"$sname"_all_somatics.csv"
	snp_vcf="vcf/snp/"$sname"_all_snps.vcf"
	snp_csv="vcf/snp/"$sname"_all_snps.csv"
			
	functional_somatic_vcf="vcf/somatic/functional/"$sname"_functional_somatics.vcf"
	functional_somatic_csv="vcf/somatic/functional/"$sname"_functional_somatics.csv"
	non_functional_somatic_vcf="vcf/somatic/non_functional/"$sname"_non_functional_somatics.vcf"
	non_functional_somatic_csv="vcf/somatic/non_functional/"$sname"_non_functional_somatics.csv"
	
	echo `date` " [$$] Creating variant report files: " $i
	# filter out homopolymer sites:
	#java -jar $SnpSift filter -f $i "(FILTER != 'HOMOPOLYMER')" > $tmp_vcf
	
	### 2 ### separate SNPs from Somatic Variants using Common, SNP, and ExAC and gnomAD Allele Frequencies:
	# NB: for some weird reason 'na SNP' returns fields with SNP=true. So use SNP='.' instead of 'na SNP' in the filter
	# NB: So check for each variable how to extract empty fields: with [variable]='.' OR [na variable]
	java -jar $SnpSift filter -f $i \
	"((na COMMON) | (COMMON=0)) & ((SNP = "false") | (SNP = '.')) & ((na PON_COUNT) | (PON_COUNT<6))" > $somatic_vcf
#	"((na COMMON) | (COMMON=0)) & ((SNP = "false") | (SNP = '.')) & ((AF_EXAC < 0.001) | (na AF_EXAC)) & ((gnomAD_AF < 0.001) | (na gnomAD_AF))" > $somatic_vcf

	java -jar $SnpSift filter -f $i -n \
	"((na COMMON) | (COMMON=0)) & ((SNP = "false") | (SNP = '.')) & ((na PON_COUNT) | (PON_COUNT<6))" > $snp_vcf
#	"((na COMMON) | (COMMON=0)) & ((SNP = "false") | (SNP = '.')) & ((AF_EXAC < 0.001) | (na AF_EXAC)) & ((gnomAD_AF < 0.001) | (na gnomAD_AF))" > $snp_vcf

	### 2 ### separate on impact: only for tumor samples; not for the PON control samples
	# eventueel TODO: CLINSIG impact: benign/ Pathogenic
	java -jar $SnpSift filter -f $somatic_vcf \
	"(ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE')" > $functional_somatic_vcf

	java -jar $SnpSift filter -f $somatic_vcf -n \
	"(ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE')" > $non_functional_somatic_vcf

	# What you may really want are lines where ANY effect to be missense_variant:
	# ANN[*].EFFECT has 'missense_variant'

####################################################################################
	### 3 ### Extract desired columns from the vcf files (see Annotation_v1 script):

	java -jar $SnpSift extractFields -e "." $somatic_vcf \
CHROM POS "ANN[0].GENE" REF ALT "GEN[*].DP" "GEN[*].FREQ" "GEN[*].ADF" "GEN[*].ADR" \
"ANN[0].IMPACT" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
ID "COMMON" "PON_COUNT" "RS" "CAF" "LOF" "NMD" "MUT" \
"CLNSIG" "ORIGIN" "SNP" "AF_EXAC" "AF_TGP" "gnomAD_AF" \
"COSM.ID" "FATHMM" "MUT.ST"> $somatic_csv

	java -jar $SnpSift extractFields -e "." $snp_vcf \
CHROM POS "ANN[0].GENE" REF ALT "GEN[*].DP" "GEN[*].FREQ" "GEN[*].ADF" "GEN[*].ADR" \
"ANN[0].IMPACT" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
ID "COMMON" "PON_COUNT" "RS" "CAF" "LOF" "NMD" "MUT" \
"CLNSIG" "ORIGIN" "SNP" "AF_EXAC" "AF_TGP" "gnomAD_AF" \
"COSM.ID" "FATHMM" "MUT.ST"> $snp_csv

	java -jar $SnpSift extractFields -e "." $functional_somatic_vcf \
CHROM POS "ANN[0].GENE" REF ALT "GEN[*].DP" "GEN[*].FREQ" "GEN[*].ADF" "GEN[*].ADR" \
"ANN[0].IMPACT" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
ID "COMMON" "PON_COUNT" "RS" "CAF" "LOF" "NMD" "MUT" \
"CLNSIG" "ORIGIN" "SNP" "AF_EXAC" "AF_TGP" "gnomAD_AF" \
"COSM.ID" "FATHMM" "MUT.ST"> $functional_somatic_csv

	java -jar $SnpSift extractFields -e "." $non_functional_somatic_vcf \
CHROM POS "ANN[0].GENE" REF ALT "GEN[*].DP" "GEN[*].FREQ" "GEN[*].ADF" "GEN[*].ADR" \
"ANN[0].IMPACT" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
ID "COMMON" "PON_COUNT" "RS" "CAF" "LOF" "NMD" "MUT" \
"CLNSIG" "ORIGIN" "SNP" "AF_EXAC" "AF_TGP" "gnomAD_AF" \
"COSM.ID" "FATHMM" "MUT.ST"> $non_functional_somatic_csv
	
	# remove tmp files
	#rm $tmp_vcf $vcf_tmp2

done



