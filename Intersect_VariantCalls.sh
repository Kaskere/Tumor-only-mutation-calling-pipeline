#!/bin/bash
######################
# Consensus calling
# create intersections from vcf files of different callers
# in this case: LoFreq and VarScan

# assumed directory structure:
# Consensus calling
	# ConsensusCalling/LoFreq/vcf/not_blacklisted	(input)
	# ConsensusCalling/VarScan/vcf/not_blacklisted	(input)
	# ConsensusCalling/vcf/intersect (output)
# So, create a directory ConsensusCalling, parallel to LoFreq and VarScan dir to perform analysis 
# wd: ConsensusCalling/

# Manually copy desired vcf files
# TODO: 2 options: 
# perform both LoFreq and VarScan from scratch
# copy desired vcf files to directory for Consensus calling

######################
# dependencies

# For use with bcftools, all vcf files need to be bgzipped and tabix indexed

#bgzip="/net/nfs/PAT/home/matias/tools/tabix-0.2.6/bgzip"
#tabix="/net/nfs/PAT/home/matias/tools/tabix-0.2.6/tabix"

#for i in ./*.vcf
#do 
#bn=`basename $i`
#sname=${bn/.vcf/}
#bgzipped=$sname".vcf.gz"
#$bgzip -c $i > $bgzipped
#$tabix -s1 -b2 -e2 $bgzipped
#done

######################

# create output dirs
if [ ! -d "vcf" ]; then mkdir "vcf" ; fi
if [ ! -d "vcf/intersect" ]; then mkdir "vcf/intersect" ; fi
if [ ! -d "vcf/outersect" ]; then mkdir "vcf/outersect" ; fi


for i in LoFreq/vcf/not_blacklisted/*_not_blacklisted.vcf.gz
do 
	
	# in- and output files:
	bn=`basename $i`
	sname=${bn/_not_blacklisted.vcf.gz/}

	#LoFreq_vcf == $i
	VarScan_vcf="VarScan/vcf/not_blacklisted/"$bn
	intersect_vcf="vcf/intersect/"$sname"_intersect.vcf"
	outersect="vcf/outersect/"$sname"_outersect"

	# create intersect: keep records in LoFreq vcf, that are also in VarScan vcf
    bcftools isec -n=2 -w1 $i $VarScan_vcf -o $intersect_vcf -O v
	# create complement: print list of all sites called by one or all tools (for testing) 
	bcftools isec $i $VarScan_vcf -c all -n +0 -o $outersect -O v

done


#	i='LoFreq/vcf/not_blacklisted/TAM15-34226FF_S21_L004_not_blacklisted.vcf.gz'
#	i='LoFreq/vcf/not_blacklisted/TVU12-11384FF_S28_L004_not_blacklisted.vcf.gz'















