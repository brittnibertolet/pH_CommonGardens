#### Indicator species analysis

rm(list=ls())

# Read in data and cleaning ####
library(vegan)
library(indicspecies)
source("~/OneDrive - nd.edu/R-functions/setTheme_BB.R")
source("code/indicatorFunctions.R")
library(tidyverse)

# Read in OTU data
otu=read.table("16S/otutab_raw.txt", sep="\t", header=T)
rownames(otu)=otu$OTU.ID
# Read in metadata
meta=read.csv("16S/metadata.csv", stringsAsFactors = F)
# Read in taxonomy
tax=read.table("16S/zotus_sintax.txt", sep="\t")

# Get rid of otus in blanks 
#get OTUs in extraction blanks 
blanks=otu[,grep("ex", colnames(otu))]
blanks=rownames(blanks)[rowSums(blanks)!=0]
#remove OTUs in blank samples
otu=otu[!(otu$OTU.ID%in%c(blanks)),]
#throw out samples that failed 
otu=otu[,colnames(otu)!="pH20_NG_N3_S41_L001"]
rm(blanks)

# Get only experimental observations
otuE=otu[,colnames(otu)%in%meta$sampleID[meta$class%in%c("experiment")]]
# rarify otuE table 
min(colSums(otuE))
otuE=rrarefy(t(otuE), sample=25523)

#Look at only methanogens
methanoOTUs=as.character(tax[grepl("Methano", tax$V2),]$V1)
otuE.M=otuE[,colnames(otuE)%in%methanoOTUs]
#throw out methanogen otus from otuE
otuE.B=otuE[,!(colnames(otuE)%in%methanoOTUs)]
# Transpose relative abundances to 1 
otuE.M.ra=decostand(otuE.M, method="total", MARGIN = 1) 
otuE.B.ra=decostand(otuE.B, method="total", MARGIN = 1) 
rowSums(otuE.M.ra)
rowSums(otuE.M.ra)

# Create pH treatment group vector 
# 1 = acidic, 2 = near-neutral
pH.groups=c(rep(1, 4), rep(2, 4), #CB
            rep(1, 4), rep(2, 4), #MI
            rep(1, 4), rep(2, 3), #NG
            rep(1, 4), rep(2, 4), #PR
            rep(1, 4), rep(2, 4), #RP
            rep(1, 4), rep(2, 4), #TF
            rep(1, 4), rep(2, 4)) #WA

# Create inoculum  group vector 
# 1 = acidic source, 2 = near-neutral source, 3 = regional pool
inoc.groups=c(rep(1, 8), #CB
              rep(1, 8), #MI
              rep(1, 7), #NG
              rep(2, 8), #PR
              rep(3, 8), #RP
              rep(2, 8), #TF
              rep(2, 8)) #WA

# Indicator species analysis across pH treatments
set.seed(13)
indvalM.pH = multipatt(otuE.M.ra, pH.groups,
                       control = how(nperm=1000))
indvalM.pH.out = extractIndicatorZOTUs.pH(data=indvalM.pH$sign)
indvalB.pH = multipatt(otuE.B.ra, pH.groups,
                       control = how(nperm=1000))
indvalB.pH.out = extractIndicatorZOTUs.pH(data=indvalB.pH$sign)

# Indicator species analysis across source inoculum
indvalM.inoc = multipatt(otuE.M.ra, inoc.groups,
                         control = how(nperm=1000))
indvalM.inoc.out = extractIndicatorZOTUs.inoc(data=indvalM.inoc$sign)
indvalB.inoc = multipatt(otuE.B.ra, inoc.groups,
                         control = how(nperm=1000))
indvalB.inoc.out = extractIndicatorZOTUs.inoc(data=indvalB.inoc$sign)

# How many unique methanogen ZOTUs associated with the treatments?
length(unique(c(indvalM.pH.out$ZOTUs, indvalM.inoc.out$ZOTUs)))
# Across pH, number associated with each pH treatment
nrow(indvalM.pH.out[indvalM.pH.out$pH=="acidic",])
nrow(indvalM.pH.out[indvalM.pH.out$pH=="near-neutral",])
# Across inoculum source, number associated with each pH type
nrow(indvalM.inoc.out[indvalM.inoc.out$pH=="acidic",])
nrow(indvalM.inoc.out[indvalM.inoc.out$pH=="near-neutral",])
# How many unique non-methanogen ZOTUs associated with the treatments?
length(unique(c(indvalB.pH.out$ZOTUs, indvalB.inoc.out$ZOTUs)))
# Across pH, number associated with each pH treatment
nrow(indvalM.pH.out[indvalB.pH.out$pH=="acidic",])
nrow(indvalM.pH.out[indvalB.pH.out$pH=="near-neutral",])
# Across inoculum source, number associated with each pH type
nrow(indvalM.inoc.out[indvalB.inoc.out$source=="acidic",])
nrow(indvalM.inoc.out[indvalB.inoc.out$source=="near-neutral",])


# Generate heat maps
plot_grid(genHeatmap.MCC(ind=indvalM.pH.out, label="pH treatment"),
          genHeatmap.MCC(ind=indvalM.inoc.out, label="inoculum source"),
          ncol=2, labels=c("A","B"))
ggsave("figures/heatMCC.pdf", height=4, width=6)

plot_grid(genHeatmap.nonMCC(ind=indvalB.pH.out, label="pH treatment"),
          genHeatmap.nonMCC(ind=indvalB.inoc.out, label="inoculum source"),
          ncol=2, labels=c("A","B"))
ggsave("figures/heatNonMCC.pdf", height=8, width=6.5)

write.csv(indvalM.pH.out$ZOTUs, file = "trees/indicator.M.pH.csv", quote = F, row.names = F)
write.csv(indvalB.pH.out$ZOTUs, file = "trees/indicator.B.pH.csv", quote = F, row.names = F)
write.csv(indvalM.inoc.out$ZOTUs, file = "trees/indicator.M.inoc.csv", quote = F, row.names = F)
write.csv(indvalB.inoc.out$ZOTUs, file = "trees/indicator.B.inoc.csv", quote = F, row.names = F)

write.csv(indvalM.pH.out, file = "trees/indicator.M.pH2.csv", quote = F, row.names = F)
write.csv(indvalB.pH.out, file = "trees/indicator.B.pH2.csv", quote = F, row.names = F)
write.csv(indvalM.inoc.out, file = "trees/indicator.M.inoc2.csv", quote = F, row.names = F)
write.csv(indvalB.inoc.out, file = "trees/indicator.B.inoc2.csv", quote = F, row.names = F)




