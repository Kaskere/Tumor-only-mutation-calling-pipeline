##############################################################################################################
##############################################################################################################
# Date 19 june 2018
# Create Mutation profiles for PELs, 1 august

# Mutation frequency bar plots:

##############################################################################################################
##############################################################################################################
### SET gene_blacklist ###

# genes described in Lawrence et al, Nature, 2013, as passenger mutations:
# gene_blacklist <- c('MUC4', 'MUC16', 'PCLO', 'AR', 'TTN', 'LRP1B', 'CSMD3')

# other putative passenger genes
#'TCF4', 'KLHL14', 'SYNE1', 'DCHS1', 'DST', 'CELSR2'
#'HIST1H2AL', 'KDM6B', 'KLHL14', 'ARID3A')

# putative oncogenic effect:
# ZFHX3, ZFHX4, IRF4, FOXO3, TSC2
# proved oncogenic effect:
# ARID1A, PIK3CA

# read results from MutSigCV on AllSomatic PEL data:
# wd: /net/nfs/PAT/analysis/MPS-339/PELLL_PE125/MutSigCV/AllSomatic
#out <- read.csv('/net/nfs/PAT/analysis/MPS-339/PELLL_PE125/MutSigCV/AllSomatic/PELn8vaf10.sig_genes.csv',  sep=',', header=T, stringsAsFactors=F, row.names=1)
#sig.genes.p0.05 <- out$gene[which(out$p<0.05)]
#pel8_vaf10_mutsig_p0.05 <- sig.genes.p0.05[-1]

# chapuy_mutsig <- read.csv('/net/nfs/PAT/analysis/MPS-339/PELLL_PE125/LoFreq/significantGenes_Chapuy_MutSig2CV.csv', sep='\t', header=T, stringsAsFactors=F)

########################################################
# Mutation dataframe containing all somatic mutations  #
########################################################
# wd: /net/nfs/PAT/analysis/MPS-339/PELLL_PE125/ConsensusCalls

# path to vcf files:
# /net/nfs/PAT/analysis/MPS-339/PELLL_PE125/ConsensusCalls/vcf/somatic/functional
vcf_path <- getwd()

# index of all csv files
mutfileIx <- grep('.csv', list.files(vcf_path) ) 
length(mutfileIx) # == 8

# number of somatic mutations per sample
All_nr_muts_sample <- c()

All_Functional_Somatic_DF <- c()

# for each somatic.csv file
for (j in 1:length(mutfileIx)){

	# Lees csv file in en sample name:
	mutfileName <- list.files(vcf_path)[mutfileIx][j]
	mutFile <- read.csv(mutfileName, sep='\t', header=T, stringsAsFactors=F)
	sampleName <- strsplit(mutfileName, '_functional_somatics.csv', fixed=T)[[1]]
	
	print(mutfileName)
	print(sampleName)
	print(j)

	# if sample has no mutation, add the sampleName + 'NA' to the mutFile
	if(nrow(mutFile)==0){
		mutFile[1,] <- rep(NA, ncol(mutFile))
	}

	# Add sampleName to each row in csv file
	Sample <- rep(sampleName, nrow(mutFile))
	mutFile <- cbind(mutFile, as.character(Sample))
	
	# Count number of mutations in sample and add to vector
	muts_in_sample <- nrow(mutFile)
	All_nr_muts_sample <- c(All_nr_muts_sample, muts_in_sample)
	
	# Add to MasterDF
	All_Functional_Somatic_DF <- rbind(All_Functional_Somatic_DF, mutFile)
}	


# Reorder data frame:
All_Functional_Somatic_DF <- All_Functional_Somatic_DF[, c(33,1:32)]
colnames(All_Functional_Somatic_DF) <- c('PA.id','CHROMO','POS','GENE','REF','ALT','DEPTH', 'VAF', 'Alt_fwd', 'Alt_rev', 'Strand bias', 'HRUN', 'IMPACT','EFFECT', 'NTchange', 'AAchange', 
'dbSNP.id','COMMON', 'PON_COUNT', 'RS', 'CAF', 'LOF', 'NMD', 'MUT', 'CLNSIG', 'ORIGIN', 'SNP', 'AF_ExAC', 'AF_TGP', 'AF_gnomAD', 'COSM.ID', 'FATHMM', 'MUT.ST')

row.names(All_Functional_Somatic_DF) <- 1:nrow(All_Functional_Somatic_DF)

#All_Functional_Somatic_DF <- AllSomatic_DF
#dir.create('data/MutationTables')

#write.table(All_Functional_Somatic_DF, file='/net/nfs/PAT/analysis/MPS-339/PELLL_PE125/LoFreq/vcf/somatic/All_Filt_High_impact_SomaticSNVs.txt', sep='\t', quote=F)
write.table(All_Functional_Somatic_DF, file="/net/nfs/PAT/analysis/MPS-339/PELLL_PE125/ConsensusCalls/vcf/somatic/All_Functional_Somatic_DF.txt", sep='\t', quote=F)
#write.table(All_Functional_Somatic_DF, file='/net/nfs/PAT/analysis/MPS-339/PELLL_PE125/LoFreq/vcf/somatic/PEL8_MutSigCVp005_AllSomaticSNVs.txt', sep='\t', quote=F)



########################################################
# 2 # create tables: MutPrevalencePerPatient.txt, MutNumbersPerPatient.txt, GeneFreqs.txt
########################################################
#All_Functional_Somatic_DF <- read.table('/net/nfs/PAT/analysis/MPS-339/PELLL_PE125/LoFreq/vcf/somatic/LawrenceGenesFilter_AllSomaticSNVs.txt', sep='\t', stringsAsFactors=F)

allSNVs <- All_Functional_Somatic_DF
#allSNVs$AF_gnomAD <- as.numeric(allSNVs$AF_gnomAD)
#allSNVs <- allSNVs[-which(allSNVs$AF_gnomAD>1.0e-04),]
#allSNVs <- allSNVs[-which(allSNVs$CLNSIG=='Benign/Likely_benign'),]
#allSNVs <- allSNVs[-which(allSNVs$CLNSIG=='Benign'),]

setwd('/net/nfs/PAT/analysis/MPS-339/PELLL_PE125/ConsensusCalls/vcf/somatic/functional')


# CODE HIERONDER HOEFT NIET AANGEPAST,
# SLECHTS HIERBOVEN allSNVs VERVANGEN MET MUTATION TABLE PER GROEP
# EN ONDERAAN ALLEEN OPLETTEN MET WEGSCHRIJVEN (heb zelf later met de hand in bash de folders gemaakt)

all.genes <- unique(as.character(allSNVs$GENE))
all.cases <- unique(as.character(allSNVs$PA.id))

# fill matrix with all mutations
MutNumbers <- matrix(data = NA, ncol = length(all.genes), nrow = length(all.cases))
colnames(MutNumbers) <- all.genes
rownames(MutNumbers) <- all.cases
MutPrevalence <- MutNumbers

# for each patient
for (i in 1:nrow(MutNumbers)){

	#print(i)
	sname <- rownames(MutNumbers)[i]
	# which rows in mutation-df correspond to this sname
	row.nr <- which(sname==allSNVs$PA.id)

	# what are the genes that are mutated for this sample (and not NA):
	ix <- !is.na(unique(allSNVs$GENE[row.nr]))
	mutgenes <- as.character(unique(allSNVs$GENE[row.nr])[ix])

	# # Number of mutations
	# which column/gene to fill:
	col.nr <- match( mutgenes, colnames(MutNumbers) )
	# if no mutation present: mark all genes for this sample with '0'
	if(length(col.nr)==0) {
		MutNumbers[i, ] <- 0
		MutPrevalence[i, ] <- 0
	} else{
	# voor elk gemuteerd gen, tel het aantal mutaties in dat gen
	mutGeneIx <- unique(col.nr[which(!is.na(col.nr))])
	for (z in 1:length(mutGeneIx)){
		count <- length(which(allSNVs$GENE[row.nr]==mutgenes[z]))
		type <- 
		#count <- length(which(col.nr==mutGeneIx[z]))
		MutNumbers[i, mutGeneIx[z]] <- count
		MutPrevalence[i, mutGeneIx[z]] <- 1
		}
	# voor alle andere genen/kolommen (geeen mutatie) <- 0
	noMut <- is.na(MutNumbers[i,])
	MutNumbers[i, noMut] <- 0
	MutPrevalence[i, noMut] <- 0
	}
}

### save table: number of mutations per gene per patient
#write.table(MutNumbers, file='PEL.MutNumbersPerPatient.txt', sep='\t', quote=F)
#write.table(MutNumbers, file='LawrenceGenesFilter_MutNumbersPerPatient.txt', sep='\t', quote=F)
write.table(MutNumbers, file='PEL_N8_MutNumbersPerPatient.txt', sep='\t', quote=F)

### save table: mutations present per patient yes/no (1/0)
#write.table(MutPrevalence, file='PEL.MutPrevalencePerPatient.txt', sep='\t', quote=F)
#write.table(MutPrevalence, file='LawrenceGenesFilter_MutPrevalencePerPatient.txt', sep='\t', quote=F)
write.table(MutPrevalence, file='PEL_N8_MutPrevalencePerPatient.txt', sep='\t', quote=F)

## Now for each GENE, get the
# min, max, and mean number of mutations within 1 patient
# and frequency in patients

stats <- c('min', 'max', 'mean', 'freq')
MutNumbStats <- matrix(data = NA, ncol = length(all.genes), nrow = length(stats))
colnames(MutNumbStats) <- all.genes
rownames(MutNumbStats) <- stats

for (i in 1:length(all.genes)){

	MutNumbStats[1,i] <- min(MutNumbers[,i])
	MutNumbStats[2,i] <- max(MutNumbers[,i])
	MutNumbStats[3,i] <- mean(MutNumbers[,i])
	MutNumbStats[4,i] <- 100*(sum(MutPrevalence[,i])/nrow(MutNumbers))

}

MutNumbStats <- t(MutNumbStats)
MutNumbStats <- MutNumbStats[ order(-MutNumbStats[,2]) ,]

### save data
#write.table(MutNumbStats, file='PEL.aSHMstats.txt', sep='\t', quote=F)
write.table(MutNumbStats, file='PEL_N8_aSHMstats.txt', sep='\t', quote=F)
#write.table(MutNumbStats, file='PEL8_MutSigCVp005_aSHMstats.txt', sep='\t', quote=F)

GENEFREQS <- MutNumbStats[ order(-MutNumbStats[,4]) ,]
GENEFREQS[,3] <- round(GENEFREQS[,3], digits=2)
GENEFREQS[,4] <- round(GENEFREQS[,4], digits=2)


### save data
#write.table(GENEFREQS, file='PEL.GeneFreqs.txt', sep='\t', quote=F)
write.table(GENEFREQS, file='PEL_N8_GeneFreqs.txt', sep='\t', quote=F)
#write.table(GENEFREQS, file='PEL8_MutSigCVp005_GeneFreqs.txt', sep='\t', quote=F)

#######################################################################
# Create list of all unique AAchanges and counts
# chek for HOTSPOTS, but also look for artifacts
#######################################################################
nrow(allSNVs) # 247

# Check alle AA changes en hoevaak ze voorkomen in de dataset:
allAAs <- unique(allSNVs$AAchange) # length(allAAs) = 229
allAAcounts <- c()

for (i in 1:length(unique(allSNVs$AAchange))){

	count <- length(grep(unique(allSNVs$AAchange)[i], allSNVs$AAchange, fixed=T))
	allAAcounts <- c(allAAcounts, count)
}
abundanceAAs <- cbind(allAAs, allAAcounts)
test <- as.data.frame(abundanceAAs, stringsAsFactors=F)
test$allAAcounts <- as.numeric(test$allAAcounts)
commonAAs <- test[order(test$allAAcounts, decreasing=T), ]

allAAs_NTs <- unique( allSNVs[ ,c(2:6,15:16) ] ) # 408

all_collapsed_mutations <- c()
# collapse all entries
for (i in 1:nrow(allSNVs)){	
	collapsed_mutation <- paste( allSNVs[i ,c(2:6,15:16) ], collapse=";")
	all_collapsed_mutations <- rbind(all_collapsed_mutations, collapsed_mutation)
}

all_NTc_counts <- c()
for (i in 1:nrow(allAAs_NTs)){

	# collapsed pattern to match
	pattern <- paste(allAAs_NTs[i,], collapse=";")
	count <- length(which(pattern==all_collapsed_mutations))
	all_NTc_counts <- c(all_NTc_counts, count)
}

test <- cbind(allAAs_NTs, all_NTc_counts)
final <- test[order(test$all_NTc_counts, decreasing=T), ]
#dots <- which(final$AAchange=='.')
#almost_there <- final[-dots, ]
#NAs <- which(is.na(almost_there$AAchange))
#FINAL <- almost_there[-NAs,]

colnames(final)[8] <- 'COUNTS'

FREQ <- 100*(final[,8]/length(unique(allSNVs$PA.id)))
FREQ <- round(FREQ, 0)
FINAL <- cbind(final, FREQ)

rownames(FINAL) <- c()

write.table(FINAL, file='all_unique_AAchanges_and_counts.txt', sep='\t', quote=F)

# take a look at all the different mutations in some genes
#FINAL[which(FINAL$GENE=='MYD88'),]
#FINAL[which(FINAL$GENE=='CREBBP'),]
#FINAL[which(FINAL$GENE=='SETD1B'),]
#FINAL[which(FINAL$GENE=='NOTCH2'),]





