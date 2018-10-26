#!/bin/bash
###############################################

# Date: 29 mei 2018
# input: dedupped, coordinate-sorted bam files (recalibrated also..?)
# 


# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # ????
# Set permissions for MPS-###/runfolder/ 

###############################################
# GLOBAL SETTINGS

# Provide directory with bam files. select bam files, based on extension (string: recal.bam / coordsorted.bam)

# number of threads to use: THREADS=

# mutation target bed file: TARGET_REGIONS_BED='/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/BCNHLv2_allExons.bed'

# minimal Variant Allele Frequency: minVAF=

# Mutation Callers to use: VarScan2, LoFreq, VarDict

# Combining VCFs yes/no. what overlap? 




###############################################
# GENERAL SETTINGS PER STEP 
 
# Variant Callers general settings
	# similar output - names/fields
	# SORT VCF

# Annotation general settings:

# Variant Discrimination general settings:
	# use minVAF
	# remove variants outside of target regions

###############################################

###############################################
# DEPENDENCIES & PREDEFINED VARIABLES
###############################################

# modules & variables
#java="/ccagc/lib/java/jre1.8.0_25/bin/java"
module load java
module load bwa
module load samtools

#fastqc="/net/nfs/PAT/lib/FastQC/0.11.2/fastqc"

manifest="/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/"
target_bed='/net/nfs/PAT/home/matias/data/manifests/BCNHL_Seq_v2/BCNHLv2_allExons.bed'

#ref="/ccagc/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
ref="/net/nfs/PAT/data/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"

#dbsnp="/ccagc/data/ref/gatk_bundle/2.8/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf"
#dbsnp="/net/nfs/PAT/data/ref/gatk_bundle/2.8/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf" 

############################################
# OPTIONAL STEP: 
# Select aligned reads from bam file that overlap with capture target regions 
# for quicker processing (for somatic mutation analysis only)
# downside: SNPs outside of target regions (introns) are not called/evaluated
# if not used, still after each calling tool, only mutations in regions of interest are kept.

if [ ! -e log/targetRegionFilter ]
then
	./code/targetRegionFilter.sh && touch log/targetRegionFilter
else
	echo "log/targetRegionFilter already exists: skipping targetRegionFilter" 
fi

############################################

if [ ! -e log/LoFreq_bam_indelQ ]
then
	./code/LoFreq_bam_indelQ.sh && touch log/LoFreq_bam_indelQ
else
	echo "log/LoFreq_bam_indelQ already exists: skipping LoFreq_bam_indelQ" 
fi


############################################

# VarScan: Variant Calling

if [ ! -e log/VarScan_VariantCalling ]
then
	./code/VarScan_VariantCalling.sh && touch log/VarScan_VariantCalling
else
	echo "log/VarScan_VariantCalling already exists: skipping VarScan_VariantCalling" 
fi

############################################

# VarScan: Blacklist filter

if [ ! -e log/VarScan_BlacklistFilter ]
then
	./code/VarScan_BlacklistFilter.sh && touch log/VarScan_BlacklistFilter
else
	echo "log/VarScan_BlacklistFilter already exists: skipping VarScan_BlacklistFilter" 
fi

############################################

# VarScan: annotation

if [ ! -e log/VarScan_Annotation ]
then
	./code/VarScan_Annotation.sh && touch log/VarScan_Annotation
else
	echo "log/VarScan_Annotation already exists: skipping VarScan_Annotation" 
fi

############################################

# Separate SNPs and Somatic variants from VCF file:
if [ ! -e log/VarScan_VariantDiscrimination ]
then
	./code/VarScan_VariantDiscrimination.sh && touch log/VarScan_VariantDiscrimination
else
	echo "log/VarScan_VariantDiscrimination already exists: skipping VarScan_VariantDiscrimination" 
fi

############################################
# COMBINE VCF FILES FROM DIFFERENT CALLERS








############################################





