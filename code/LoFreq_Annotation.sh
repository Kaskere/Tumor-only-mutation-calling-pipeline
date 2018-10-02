#!/bin/bash

# Variant Annotation
# NB: Use uncompressed VCF files/databases 

module load java
#snpEff="/net/nfs/PAT/home/matias/tools/snpEff_4.1/"
snpEff="/net/nfs/PAT/home/matias/tools/snpEff_4.3/"

tabix="/net/nfs/PAT/home/matias/tools/tabix-0.2.6/tabix"
bgzip="/net/nfs/PAT/home/matias/tools/tabix-0.2.6/bgzip"
ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
manifest="/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/"
#dbSnp="/net/nfs/PAT/home/matias/dbSNP/b150/common_all_20170710.vcf"
dbSnp="/net/nfs/PAT/home/matias/dbSNP/b151/All_20180423.vcf.gz"
ClinVar="/net/nfs/PAT/home/matias/dbSNP/b151/clinvar_20180701.vcf.gz"
Cosmic="/net/nfs/PAT/home/matias/data/ref/cosmic/hg19_v84_2018/CosmicCodingMuts_BCNHLv2_v84_hg19.vcf"
#PON_blacklist="/net/nfs/PAT/home/matias/data/blacklist/complete.PON_blacklist.vcf"
#PON_blacklist="/net/nfs/PAT/home/matias/data/blacklist/blacklist_PON22_LoFreq.vcf"
#PON_blacklist="/net/nfs/PAT/home/matias/data/blacklist/blacklist_PON22_LoFreq_RAW.vcf.gz"
PON_blacklist="/net/nfs/PAT/home/matias/data/blacklist/blacklist_PON24_LoFreq_RAW.vcf"
#gnomAD='/net/nfs/PAT/home/tjitske/dbSNP/gnomAD/2.0.2/retagged_af-only-gnomad.raw.sites.b37.vcf.gz'
gnomAD='/net/nfs/PAT/home/tjitske/dbSNP/gnomAD/2.0.2/retagged_gnomad.exomes.r2.0.2-AF.vcf.gz'

#mkdir vcf/annotated
if [ ! -d "vcf/annotated" ]; then mkdir "vcf/annotated" ; fi	
if [ ! -d "vcf/annotated/stats" ]; then mkdir "vcf/annotated/stats" ; fi	

# Variant Annotation  
for i in vcf/raw/*_lofreq.vcf
do
	# in- and output files:
	bn=`basename $i`
	sname=${bn/_lofreq.vcf/}
	snpEff_stats="vcf/annotated/"$sname"_stats.html"
	vcf_tmp="vcf/annotated/"$sname"_tmp.vcf"
	vcf_tmp2="vcf/annotated/"$sname"_tmp2.vcf"
	vcf_tmp3="vcf/annotated/"$sname"_tmp3.vcf"

	vcf_out="vcf/annotated/"$sname"_annot.vcf"
	vcf_out_gz="vcf/annotated/"$sname"_annot.vcf.gz"
	vcf_csv="vcf/annotated/"$sname"_annot.csv"

	###### dbSNP annotation, Clinvar annotation, Effect prediction
	#if [ ! -e $vcf_tmp ]; then do
	java -Xmx4g -jar $snpEff/SnpSift.jar annotate $dbSnp $i \
	| java -Xmx4g -jar $snpEff/SnpSift.jar annotate $ClinVar - \
	| java -Xmx4g -jar $snpEff/snpEff.jar eff -filterInterval $manifest/BCNHLv2_primary_coord.bed -v -canon -strict -stats $snpEff_stats hg19 - > $vcf_tmp
	#; fi

	# | java -Xmx4g -jar $snpEff/SnpSift.jar annotate -id $Cosmic - \

	### COSMIC, gnomAD and Panel of Normal annotation
	java -Xmx4g -jar $snpEff/SnpSift.jar annotate -v $Cosmic $vcf_tmp \
	| java -Xmx4g -jar $snpEff/SnpSift.jar annotate -v -info 'gnomAD_AF' $gnomAD - \
	| java -Xmx4g -jar $snpEff/SnpSift.jar annotate -v -noAlt -info 'PON_BLACKLIST' $PON_blacklist - > $vcf_out


	### Create csv files - extract fields of interest
java -jar $snpEff/SnpSift.jar extractFields -e "." $vcf_out \
CHROM POS "ANN[0].GENE" REF ALT DP AF \
"ANN[0].IMPACT" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
ID "COMMON" "RS" "CAF" "LOF" "NMD" "MUT" \
"CLNSIG" "ORIGIN" "SNP" "AF_EXAC" "AF_TGP" "gnomAD_AF" \
"COSM.ID" "FATHMM" "MUT.ST" "PON_BLACKLIST"> $vcf_csv

	# remove tmp files
	rm $vcf_tmp
	#$vcf_tmp2 $vcf_tmp3

	# mv stats files
	mv vcf/annotated/*_stats.genes.txt vcf/annotated/stats/
	mv vcf/annotated/*_stats.html vcf/annotated/stats/

	# index the remaining vcf files:
	$bgzip -c $vcf_out > $vcf_out_gz
	$tabix -p vcf $vcf_out_gz
	
	
done
# EOF


### COSMIC annotation
#java -jar $snpEff/SnpSift.jar annotate -v $Cosmic $vcf_tmp > $vcf_tmp2
### gnomAD annotation
#java -jar $snpEff/SnpSift.jar annotate -v -info 'gnomAD_AF' $gnomAD $vcf_tmp2 > $vcf_tmp3
### Panel of Normal annotation (only add the PON_BLACKLIST tag; blacklisted y/n): NB position based (-noAlt).				
#java -jar $snpEff/SnpSift.jar annotate -v -noAlt -info 'PON_BLACKLIST' $PON_blacklist $vcf_tmp3 > $vcf_out


#i='/net/nfs/PAT/analysis/MPS-310/SNVanalysis/LoFreq_dev3/vcf/somatic/TAM15-34226_S22_L004_somatics.vcf'

#tabix="/net/nfs/PAT/home/matias/tools/tabix-0.2.6/tabix"
#bgzip="/net/nfs/PAT/home/matias/tools/tabix-0.2.6/bgzip"

#$bgzip -c $i > TAM15.vcf.gz
#$tabix -s1 -b2 -e2 TAM15.vcf.gz

#bcftools view -O v TAM15.vcf.gz \
#	| sed -e 's/\([;=[:space:]]\)AF\([,;=[:space:]]\)/\VAF\2/' > TAM15_retagged.vcf.gz
#  | bcftools view -O z -o variants.re-tagged.vcf.gz



