# *Onthophagus taurus* population transcriptomics study

This repository holds scripts for the data processing and analysis for a common garden transcriptomics experiment of three *Onthophagus taurus* populations. 

One native (Italy, IT) and two exotic (Western Australia, WA and Eastern US, North Carolina, NC) were kept in common garden conditions for three generations. 
RNA was extracted and sequenced for each of three individuals from each of four developmental stages and both sexes, resulting in 72 individuals sampled.

## Data availability

Sequence data: Raw reads (1.35 billion 2 x 100 bp paired-end Illumina) will be available at NCBI SRA: BioProject ID: TBD

## Programs used

- Trimmomatic v. 0.33
- R version 3.3.2
    + DESeq2 1.14.1
- VCFtools version 0.1.15
- PGDSpider 2.0.9.0
- BayeScan2.1 

Program versions can also be found in the scripts and/or in the associated manuscript (submitted/bioRxiv)

## Scripts

	- `Data_processing_commands.md`
		- Clean and map reads
		- Alignment, variant calling, filtering
		- Population genetic calculations

	- `DGEforPopxSex_eachStage.R`
		- Differential gene expression analysis
