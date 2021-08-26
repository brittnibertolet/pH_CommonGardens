# pH_CommonGardens

Data and code associated with "Microbial community composition, and not pH, influences lake sediment function" by Brittni L. Bertolet, Sydney I. Louden, and Stuart E. Jones. 


Data:
- rawreads.zip: Zipped directory of fastq files containing raw sequence reads for each sample 
- metadata.csv: Metadata describing samples. 
- zotutab_raw.txt: Abundance table of Zero-radius Operational Taxonomic Units (ZOTUs)
- zotus.repSequences.fa: Representative sequences for ZOTUs
- zotus_sintax.tax: Taxonomic associations with ZOTUs
- CH4CO2rates.csv: Rates of carbon dioxide and methane production associated with experimental microcosms

Code:
- community-Analyses_2021-02-17.R: R code for all initial community analyses
- indicatorFunctions.R: R functions for indicator analyses
- indicatorSpecies_2021-06-08.R: R code for conducting indicator species analysis
- phylogeneticSignal.R: R code to determine phylogenetic signal across indicator species
- usearch-CommunityAnalyses.txt: Bash code for determining ZOTUs

