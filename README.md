# LoFreq-TumorOnly
Workflow for SNV calling on tumor-only targeted NGS data

Currently, there are 2 preprocessing steps on the bam files
(but eventually these could be moved to an 'post-alignment' processing workflow) 
- target region filter
- LoFreq indelquality calculation

Includes different steps:
1. LoFreq variant calling
- followed by 'lofreq filter' based on NGS quality metrics
2. Apply LoFreq Panel Of Normal (PON) Blacklist
- to remove all artificial/germline calls from the PON dataset
(todo: determine what cutoff to use: filter out all calls present in at least 2 PON-cases)
3. Variant Annotation
4. Variant Discrimination
- distinguish snp/somatic
- optionally: from all somatic calls: distinguish driver/passenger (using MutSig2CV)
- distinguish functional/non-functional somatic
- optionally: from all functional somatic SNVs, select only the significant driver genes, based on list of significant driver genes from MutSig2CV-output


