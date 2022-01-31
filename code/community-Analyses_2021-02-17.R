rm(list=ls())

# Read in data and cleaning ####
library(vegan)
library(effectsize)
source("~/OneDrive - nd.edu/R-functions/setTheme_BB.R")

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

# Look at community composition of both experimental and inoculum sources 
otuI=otu[,colnames(otu)%in%meta$sampleID[meta$class%in%c("initial", "experiment")]]
# rarify otuI table 
min(colSums(otuI))
otuI=rrarefy(t(otuI), sample=25523)
#otuI=rrarefy(t(otuI), sample=29031)

rowSums(otuI)

#Look at only methanogens
methanoOTUs=as.character(tax[grepl("Methano", tax$V2),]$V1)
motuI=otuI[,colnames(otuI)%in%methanoOTUs]
#throw out methanogen otus from otuI
otuI=otuI[,!(colnames(otuI)%in%methanoOTUs)]
# Transpore relative abundances to 1 
otuI.ra=decostand(otuI, method="total", MARGIN = 1) 
MotuI.ra=decostand(motuI, method="total", MARGIN = 1) 
rowSums(otuI.ra)
rowSums(MotuI.ra)

#Generate distance matrix
distM=vegdist(MotuI.ra, method="bray")
distB=vegdist(otuI.ra, method="bray")

#Generate principle coordinate scaling with eigenvalues
pcoaM=cmdscale(distM, eig=T)
pcoaB=cmdscale(distB, eig=T)
#Create new dataframe 
pcoaI.DF=data.frame(sampleID=rownames(pcoaM$points),
                    Maxis1=pcoaM$points[,1],
                    Maxis2=pcoaM$points[,2],
                    nonaxis1=pcoaB$points[,1],
                    nonaxis2=pcoaB$points[,2])
rownames(pcoaI.DF)=NULL

#merge the two dataframes
df=merge(meta,pcoaI.DF, by="sampleID")
df$lakeType[df$lakeID%in%c("MI", "CB", "NG")]="acidic"
df$lakeType[df$lakeID%in%c("WA", "TF", "PR")]="near-neutral"
df$lakeType[df$lakeID%in%c("RP")]="regional pool"



#PERMANOVA
adonis(distB~df$class)
adonis(distM~df$class)

#plot MCC principle coordinates
df$treatment[df$class=="initial"]="initial"
df$treatment[df$treatment=="neutral"]="near-neutral"

df$treatment=factor(df$treatment, levels=c("acidic", "near-neutral", "initial"),
                    labels=c("acidic incubation", "near-neutral incubation", "initial source"))

motuI.p=ggplot(df, aes(x=Maxis1, y=Maxis2))+
  geom_point(aes(color=lakeType, shape=treatment))+
  xlab("PCoA Axis 1 [28.7%]")+
  ylab("\nPCoA Axis 2 [14.5%]")+
  ggtitle("MCC")+
  scale_color_discrete(name="Inoculum Source")+
  scale_shape_manual(name="Sample Type", values = c(16,17, 3))+
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.title=element_text(size=10), 
        legend.text= element_text(size=9),
        plot.title = element_text(size=10, hjust = 0.5))+
  #scale_x_continuous(limits = c(-0.45, 0.55), breaks = c(-0.25, 0, 0.25))+
  #scale_y_continuous(limits = c(-0.55, 0.3), breaks = c(-0.4,-0.2, 0, 0.2))+
  NULL
motuI.p
#look at percent explained
pcoaM$eig[1]/sum(pcoaM$eig)
pcoaM$eig[2]/sum(pcoaM$eig)

## Look at centroid of initial vs incubated:
#non-MCC
centers=aggregate(nonaxis1~lakeID+treatment, data=df, FUN="mean")
centers$nonaxis2=aggregate(nonaxis2~lakeID+treatment, data=df, FUN="mean")[,3]
centers$Maxis1=aggregate(Maxis1~lakeID+treatment, data=df, FUN="mean")[,3]
centers$Maxis2=aggregate(Maxis2~lakeID+treatment, data=df, FUN="mean")[,3]

rownames(centers)=paste(centers$lakeID, centers$treatment, sep="_")
dist(centers[centers$lakeID=="CB",3:4], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="NG",3:4], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="MI",3:4], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="PR",3:4], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="TF",3:4], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="WA",3:4], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="RP",3:4], method = "euclidean", diag = FALSE)

dist(centers[centers$lakeID=="CB",5:6], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="NG",5:6], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="MI",5:6], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="PR",5:6], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="TF",5:6], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="WA",5:6], method = "euclidean", diag = FALSE)
dist(centers[centers$lakeID=="RP",5:6], method = "euclidean", diag = FALSE)


homeDist=read.csv("distanceToHome.csv", stringsAsFactors = F)
t.test(distanceInitialToHome~lakeType, data = homeDist[homeDist$community=="MCC",])
t.test(distanceInitialToHome~lakeType, data = homeDist[homeDist$community=="non-MCC",])

ggplot(homeDist, aes(x=lakeType, y=distanceInitialToHome))+geom_boxplot()+
  facet_grid(~community)+
  ylab("Distance between centroids\n(initial - incubated in home environment)")
ggsave("distances.pdf", height=4, width=5)

#plot non-MCC principle coordinates
otuI.p=ggplot(df, aes(x=nonaxis1, y=nonaxis2))+
  geom_point(aes(color=lakeType, shape=treatment))+
  annotate("point", x=-0.01, y=-0.08)+
  annotate("point", x=0.11, y=-0.22)+
  xlab("PCoA Axis 1 [15.6%]")+
  ylab("\nPCoA Axis 2 [10.6%]")+
  ggtitle("non-MCC")+
  scale_color_discrete(name="Inoculum Source")+
  scale_shape_manual(name="Sample Type", values = c(16,17, 3))+
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.title=element_text(size=10), 
        legend.text= element_text(size=9),
        plot.title = element_text(size=10, hjust = 0.5))+
  #scale_x_continuous(limits = c(-0.45, 0.55), breaks = c(-0.25, 0, 0.25))+
  #scale_y_continuous(limits = c(-0.55, 0.3), breaks = c(-0.4,-0.2, 0, 0.2))+
  NULL
otuI.p
#look at percent explained
pcoa$eig[1]/sum(pcoa$eig)
pcoa$eig[2]/sum(pcoa$eig)


plot=plot_grid(motuI.p+guides(color=F, shape=F), get_legend(otuI.p),
               otuI.p+guides(color=F, shape=F), 
               labels=c("A", "", "B"), align = "hv", axis="t", nrow=2,
               rel_widths = c(1, 0.6))
ggsave("figures/Suppinitial-community.pdf", plot,height=5, width=5)


# Look at community composition in the inocula source ####
#Get only initial inocula source observations
otuI=otu[,colnames(otu)%in%meta$sampleID[meta$class=="initial"]]
# rarify otuI table 
min(colSums(otuI))
otuI=rrarefy(t(otuI), sample=39568)
rowSums(otuI)

#Look at only methanogens
methanoOTUs=as.character(tax[grepl("Methano", tax$V2),]$V1)
motuI=otuI[,colnames(otuI)%in%methanoOTUs]
# Transpore relative abundances to 1 
otuI.ra=decostand(otuI, method="total", MARGIN = 1) 
MotuI.ra=decostand(motuI, method="total", MARGIN = 1) 
rowSums(otuI.ra)
rowSums(motuI.ra)

#Generate distance matrix
distM=vegdist(MotuI.ra, method="bray")
dist=vegdist(otuI.ra, method="bray")

#Generate principle coordinate scaling with eigenvalues
pcoaM=cmdscale(distM, eig=T)
pcoa=cmdscale(dist, eig=T)

#Create new dataframe 
pcoaI.DF=data.frame(sampleID=rownames(pcoaM$points),
                   Maxis1=pcoaM$points[,1],
                   Maxis2=pcoaM$points[,2],
                   nonaxis1=pcoa$points[,1],
                   nonaxis2=pcoa$points[,2])
rownames(pcoaI.DF)=NULL

#merge the two dataframes
df=merge(meta,pcoaI.DF, by="sampleID")
df$lakeType[df$lakeID%in%c("MI", "CB", "NG")]="acidic"
df$lakeType[df$lakeID%in%c("WA", "TF", "PR")]="near-neutral"
df$lakeType[df$lakeID%in%c("RP")]="regional pool"

#plot MCC principle coordinates
motuI.p=ggplot(df, aes(x=Maxis1, y=Maxis2))+
  geom_point(aes(color=lakeType))+
  xlab("PCoA Axis 1 [47.8%]")+
  ylab("\nPCoA Axis 2 [24.0%]")+
  ggtitle("MCC")+
  scale_color_discrete(name="Inoculum Source")+
  scale_shape_discrete(name="pH Environment")+
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.title=element_text(size=10), 
        legend.text= element_text(size=9),
        plot.title = element_text(size=10, hjust = 0.5))+
  #scale_x_continuous(limits = c(-0.45, 0.55), breaks = c(-0.25, 0, 0.25))+
  #scale_y_continuous(limits = c(-0.55, 0.3), breaks = c(-0.4,-0.2, 0, 0.2))+
  NULL
motuI.p
#look at percent explained
pcoaM$eig[1]/sum(pcoaM$eig)
pcoaM$eig[2]/sum(pcoaM$eig)
#plot MCC principle coordinates
otuI.p=ggplot(df, aes(x=nonaxis1, y=nonaxis2))+
  geom_point(aes(color=lakeType))+
  xlab("PCoA Axis 1 [45.4%]")+
  ylab("\nPCoA Axis 2 [17.7%]")+
  ggtitle("non-MCC")+
  scale_color_discrete(name="Inoculum Source")+
  scale_shape_discrete(name="pH Environment")+
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.title=element_text(size=10), 
        legend.text= element_text(size=9),
        plot.title = element_text(size=10, hjust = 0.5))+
  #scale_x_continuous(limits = c(-0.45, 0.55), breaks = c(-0.25, 0, 0.25))+
  #scale_y_continuous(limits = c(-0.55, 0.3), breaks = c(-0.4,-0.2, 0, 0.2))+
  NULL
otuI.p
#look at percent explained
pcoa$eig[1]/sum(pcoa$eig)
pcoa$eig[2]/sum(pcoa$eig)

plot_grid(motuI.p+guides(color=F),
          otuI.p+guides(color=F),
          get_legend(otuI.p), 
          rel_widths = c(1,1,0.4), nrow=1)
ggsave("figures/initial-community.pdf", height=3, width=6.5)

# Look at methanogen community composition
MotuI.tax=as.data.frame(t(MotuI.ra))
MotuI.tax$V1=rownames(MotuI.tax)
MotuI.tax=merge(MotuI.tax, tax, by="V1")
MotuI.tax=MotuI.tax[,-10:-11]
MotuI.tax$V2=gsub("\\([01].[0-9]{0,4}\\)", "", MotuI.tax$V2)
MotuI.tax=separate(MotuI.tax, col=V2, sep=",", into=c("domain", "phylum", "class", "order", "family", "genus"))
MotuI.tax$order_genus=paste(MotuI.tax$order, MotuI.tax$genus, sep="_")
MotuI.tax$order_genus=gsub("([og]:)", "", MotuI.tax$order_genus)
MotuI.tax=MotuI.tax[,c(2:8,15)]
MotuI.tax=aggregate(.~order_genus, data=MotuI.tax, FUN = sum)
rownames(MotuI.tax)=MotuI.tax$order_genus
MotuI.tax$order_genus=NULL
MotuI.tax=as.data.frame(t(MotuI.tax))
MotuI.tax$sampleID=rownames(MotuI.tax)
MotuI.tax.melt=melt(MotuI.tax, id.vars = "sampleID", variable.name = "order_genus", value.name = "relAbund")
MotuI.tax.melt=merge(MotuI.tax.melt, df, by="sampleID")
MotuI.tax.melt$order_genus=as.character(MotuI.tax.melt$order_genus)

lowMs=colnames(MotuI.tax)[colSums(MotuI.tax[,-ncol(MotuI.tax)])/sum(MotuI.tax[,-ncol(MotuI.tax)])<0.01]
MotuI.tax.melt$order_genus[MotuI.tax.melt$order_genus%in%lowMs]="other"
MotuI.tax.melt$class
colors=c("#E56E94","#FAAFBE", "#C48793",
         #"#4E9258",
         "#786D5F",
         "#FFF380",
         "#2B547E","#157DEC",#"#82CAFA",
         "#FF7F50", 
         "#CB4335", 
         "#BDC3C7")
compI=ggplot(MotuI.tax.melt, aes(x=lakeID, y=relAbund))+
  geom_col(aes(fill=order_genus))+
  scale_fill_manual(name="Genus",values = colors, na.value="grey")+
  facet_grid(class~lakeType, space="free", scales = "free")+
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        axis.title.x = element_blank(),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.justification = "top", 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9))+
  ylab("Relative Abundance")
compI



# Look at methanogen community composition in exp ####
#Get only experimental observations
otuE=otu[,colnames(otu)%in%meta$sampleID[meta$class=="experiment"]]
#Throw out regional pool
otuE=otuE[,!colnames(otuE)%in%meta$sampleID[meta$lakeID=="RP"]]

# rarify otuE table 
min(colSums(otuE))
otuE=rrarefy(t(otuE), sample=25523)
rowSums(otuE)

#Look at only methanogens
methanoOTUs=as.character(tax[grepl("Methano", tax$V2),]$V1)
motuE=otuE[,colnames(otuE)%in%methanoOTUs]
# Transpore relative abundances to 1 
otuE.ra=decostand(otuE, method="total", MARGIN = 1) 
MotuE.ra=decostand(motuE, method="total", MARGIN = 1) 
rowSums(otuE.ra)
rowSums(MotuE.ra)

#Generate distance matrix
distM=vegdist(MotuE.ra, method="bray")
dist=vegdist(otuE.ra, method="bray")

#Generate principle coordinate scaling with eigenvalues
pcoaM=cmdscale(distM, eig=T)
pcoa=cmdscale(dist, eig=T)

#Create new dataframe 
pcoaE.DF=data.frame(sampleID=rownames(pcoaM$points),
                    Maxis1=pcoaM$points[,1],
                    Maxis2=pcoaM$points[,2],
                    nonaxis1=pcoa$points[,1],
                    nonaxis2=pcoa$points[,2])
rownames(pcoaE.DF)=NULL

#merge the two dataframes
dfE=merge(meta,pcoaE.DF, by="sampleID")
dfE$lakeType[dfE$lakeID%in%c("MI", "CB", "NG")]="acidic"
dfE$lakeType[dfE$lakeID%in%c("WA", "TF", "PR")]="near-neutral"
dfE$lakeType[dfE$lakeID%in%c("RP")]="regional pool"

#plot MCC principle coordinates
motuE.p=ggplot(dfE, aes(x=Maxis1, y=Maxis2))+
  geom_point(aes(color=lakeType, shape=treatment))+
  xlab("PCoA Axis 1 [32.1%]")+
  ylab("\nPCoA Axis 2 [16.0%]")+
  ggtitle("MCC")+
  scale_color_discrete(name="Inoculum Source")+
  scale_shape_discrete(name="pH Environment")+
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.title=element_text(size=10), 
        legend.text= element_text(size=9),
        plot.title = element_text(size=10, hjust = 0.5))+
  #scale_x_continuous(limits = c(-0.45, 0.55), breaks = c(-0.25, 0, 0.25))+
  #scale_y_continuous(limits = c(-0.55, 0.3), breaks = c(-0.4,-0.2, 0, 0.2))+
  NULL
motuE.p
#look at percent explained
pcoaM$eig[1]/sum(pcoaM$eig)
pcoaM$eig[2]/sum(pcoaM$eig)

# look at significance of treatments on methanogen community 
adonis(distM~treatment*lakeType,data=df) #sequential

fit1=adonis(distM~treatment,data=df)
fit2=adonis(distM~lakeType,data=df)
fit3=adonis2(distM~treatment*lakeType,data=df, by="margin")
#Treatment effect size = 0.19
F_to_eta2(f=fit1$aov.tab$F.Model[1], df=1, df_error = 50)
#Lake type effect size = 0.28
F_to_eta2(f=fit2$aov.tab$F.Model[1], df=2,df_error = 50)
#Interaction effect size = 0.11
F_to_eta2(f=fit3$F[1], df=2,df_error = 50)


#plot non-MCC principle coordinates
otuE.p=ggplot(dfE, aes(x=nonaxis1, y=nonaxis2))+
  geom_point(aes(color=lakeType, shape=treatment))+
  xlab("PCoA Axis 1 [18.3%]")+
  ylab("\nPCoA Axis 2 [12.0%]")+
  ggtitle("non-MCC")+
  scale_color_discrete(name="Inoculum Source")+
  scale_shape_discrete(name="pH Environment")+
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.title=element_text(size=10), 
        legend.text= element_text(size=9),
        plot.title = element_text(size=10, hjust = 0.5))+
  #scale_x_continuous(limits = c(-0.45, 0.55), breaks = c(-0.25, 0, 0.25))+
  #scale_y_continuous(limits = c(-0.55, 0.3), breaks = c(-0.4,-0.2, 0, 0.2))+
  NULL
otuE.p
#look at percent explained
pcoa$eig[1]/sum(pcoa$eig)
pcoa$eig[2]/sum(pcoa$eig)

# look at significance of treatments on non-methanogen community 
fit=adonis(dist~df$treatment*df$lakeType)

fit1=adonis(dist~treatment,data=df)
fit2=adonis(dist~lakeType,data=df)
fit3=adonis2(dist~treatment*lakeType,data=df, by="margin")

#Treatment effect size = 0.17
F_to_eta2(f=fit1$aov.tab$F.Model[1], df=1, df_error = 50)
#Lake type effect size = 0.19
F_to_eta2(f=fit2$aov.tab$F.Model[1], df=2,df_error = 50)
#Interaction effect size = 0.13
F_to_eta2(f=fit3$F[1], df=2,df_error = 50)


ggsave("figures/exp-community.pdf", height=2.5, width=3.5)

# Look at methanogen community composition
MotuE.tax=as.data.frame(t(MotuE.ra))
MotuE.tax$V1=rownames(MotuE.tax)
MotuE.tax=merge(MotuE.tax, tax, by="V1")
MotuE.tax$V3=NULL;MotuE.tax$V4=NULL
MotuE.tax$V2=gsub("\\([01].[0-9]{0,4}\\)", "", MotuE.tax$V2)
MotuE.tax=separate(MotuE.tax, col=V2, sep=",", into=c("domain", "phylum", "class", "order", "family", "genus"))
MotuE.tax$order_genus=paste(MotuE.tax$order, MotuE.tax$genus, sep="_")
MotuE.tax$order_genus=gsub("([og]:)", "", MotuE.tax$order_genus)
MotuE.tax$domain=NULL
MotuE.tax$phylum=NULL
MotuE.tax$class=NULL
MotuE.tax$order=NULL
MotuE.tax$family=NULL
MotuE.tax$genus=NULL
MotuE.tax$V1=NULL
MotuE.tax=aggregate(.~order_genus, data=MotuE.tax, FUN = sum)
rownames(MotuE.tax)=MotuE.tax$order_genus
MotuE.tax$order_genus=NULL
MotuE.tax=as.data.frame(t(MotuE.tax))
MotuE.tax$sampleID=rownames(MotuE.tax)
MotuE.tax.melt=melt(MotuE.tax, id.vars = "sampleID", variable.name = "order_genus", value.name = "relAbund")
MotuE.tax.melt=merge(MotuE.tax.melt, dfE, by="sampleID")
MotuE.tax.melt$order_genus=as.character(MotuE.tax.melt$order_genus)

lowMs=colnames(MotuE.tax)[colSums(MotuE.tax[,-ncol(MotuE.tax)])/sum(MotuE.tax[,-ncol(MotuE.tax)])<0.01]


MotuE.tax.melt$order_genus[MotuE.tax.melt$order_genus%in%lowMs]="other"
MotuE.tax.melt$bottle=paste(MotuE.tax.melt$lakeID,MotuE.tax.melt$replicate, sep=" ")
colorsE=c("#E56E94","#FAAFBE", #"#C48793",
          "#4E9258",
          #"#786D5F",
          "#FFF380",
          #"#2B547E",
          "#157DEC","#82CAFA",
          "#FF7F50", 
          "#CB4335", 
          "#BDC3C7")
MotuE.tax.melt$treatment[MotuE.tax.melt$treatment=="neutral"]="near-neutral"
compE=ggplot(MotuE.tax.melt, aes(x=bottle, y=relAbund))+
  geom_col(aes(fill=order_genus))+
  scale_fill_manual(name="Genus",values = colorsE, na.value="grey")+
  facet_grid(treatment~lakeType, space="free", scales = "free")+
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        axis.title.x = element_blank(),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.position = "top", 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9))+
  guides(fill = guide_legend(title.position = "top",
                             nrow=5, ncol= 2)) +
  scale_y_continuous(breaks=c(0, 0.5, 1))+
  ylab("Relative Abundance")
compE
ggsave("figures/expComp-2021-02-17.pdf", compE, height=4, width=6.5)


colorsL=c("#E56E94","#FAAFBE", "#C48793",
         "#4E9258",
         "#786D5F",
         "#FFF380",
         "#2B547E","#157DEC","#82CAFA",
         "#FF7F50", 
         "#CB4335", 
         "#BDC3C7")

legendDF=data.frame(val=1,genera=unique(c(unique(MotuE.tax.melt$order_genus), unique(MotuI.tax.melt$order_genus))))
l=ggplot(legendDF, aes(x=genera, y=val, fill=genera))+
  scale_fill_manual(name="Genus",values = colorsL)+geom_col()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        #legend.justification = "top", 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9))


tcomp=plot_grid(plot_grid(compI+guides(fill=F), compE+guides(fill=F), ncol=1, align = "hv", axis="l"),
          get_legend(l), ncol=1, rel_heights = c(1, 0.35))

ggsave("figures/methanoComp-2021-02-17.pdf", tcomp, height=7, width=6.5)



# Look at non-methanogen community composition in exp ####

# Look at methanogen community composition
otuE.tax=as.data.frame(t(otuE.ra))
otuE.tax$V1=rownames(otuE.tax)
otuE.tax=merge(otuE.tax, tax, by="V1")
otuE.tax$V3=NULL;otuE.tax$V4=NULL
otuE.tax$V2=gsub("\\([01].[0-9]{0,4}\\)", "", otuE.tax$V2)
otuE.tax=separate(otuE.tax, col=V2, sep=",", into=c("domain", "phylum", "class", "order", "family", "genus"))
otuE.tax$domain_phylum=paste(otuE.tax$domain, otuE.tax$phylum, sep="_")
otuE.tax$domain_phylum=gsub("([dp]:)", "", otuE.tax$domain_phylum)
otuE.tax$domain=NULL
otuE.tax$phylum=NULL
otuE.tax$class=NULL
otuE.tax$order=NULL
otuE.tax$family=NULL
otuE.tax$genus=NULL
otuE.tax$V1=NULL
otuE.tax=aggregate(.~domain_phylum, data=otuE.tax, FUN = sum)
rownames(otuE.tax)=otuE.tax$domain_phylum
otuE.tax$domain_phylum=NULL
otuE.tax=as.data.frame(t(otuE.tax))
otuE.tax$sampleID=rownames(otuE.tax)
otuE.tax.melt=melt(otuE.tax, id.vars = "sampleID", variable.name = "domain_phylum", value.name = "relAbund")
otuE.tax.melt=merge(otuE.tax.melt, dfE, by="sampleID")
otuE.tax.melt$domain_phylum=as.character(otuE.tax.melt$domain_phylum)

lowMs=colnames(otuE.tax)[colSums(otuE.tax[,-ncol(otuE.tax)])/sum(otuE.tax[,-ncol(otuE.tax)])<0.01]


otuE.tax.melt$domain_phylum[otuE.tax.melt$domain_phylum%in%lowMs]="other"
otuE.tax.melt$bottle=paste(otuE.tax.melt$lakeID,otuE.tax.melt$replicate, sep=" ")
colorsE=c("#E56E94","#FAAFBE", #"#C48793",
          "#4E9258",
          #"#786D5F",
          "#FFF380",
          #"#2B547E",
          "#157DEC","#82CAFA",
          "#FF7F50", 
          "#CB4335", 
          "#BDC3C7")
otuE.tax.melt$treatment[otuE.tax.melt$treatment=="neutral"]="near-neutral"
BcompE=ggplot(otuE.tax.melt, aes(x=bottle, y=relAbund))+
  geom_col(aes(fill=domain_phylum))+
  #scale_fill_manual(name="Genus",values = colorsE, na.value="grey")+
  facet_grid(treatment~lakeType, space="free", scales = "free")+
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        axis.title.x = element_blank(),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.position = "top", 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9))+
  guides(fill = guide_legend(title.position = "top",
                            ncol= 2)) +
  ylab("Relative Abundance")
BcompE
ggsave("figures/BexpComp-2021-02-17.pdf", BcompE, height=4, width=6.5)


# look at phya specifically

abundP=as.data.frame(colSums(otuE.tax[,-ncol(otuE.tax)])/sum(otuE.tax[,-ncol(otuE.tax)]))
colnames(abundP)="relAbund"
abundP$tax=rownames(abundP)
rownames(abundP)=NULL
fiveP=abundP[abundP$relAbund>0.10,]

# Get five most dominant
domP=otuE.tax[,colnames(otuE.tax)%in%fiveP$tax,]
domP$sampleID=rownames(domP)
domP=merge(dfE, domP, by="sampleID")


fit=aov(Bacteria_Acidobacteria~treatment*lakeType, data=domP[domP$lakeID!="RP",])
summary(fit)

fit=aov(Bacteria_Bacteroidetes~treatment*lakeType, data=domP[domP$lakeID!="RP",])
summary(fit)

fit=aov(Bacteria_Firmicutes~treatment*lakeType,data=domP[domP$lakeID!="RP",])
summary(fit)

fit=aov(Bacteria_Proteobacteria~treatment*lakeType, data=domP[domP$lakeID!="RP",])
summary(fit)


otuE.tax=merge(otuE.tax, dfE, by="sampleID")
ggplot(otuE.tax, aes(x=treatment, y=Bacteria_Proteobacteria))+
  geom_boxplot(aes(fill=lakeType))
fit=aov(Bacteria_Proteobacteria~treatment*lakeType, data=otuE.tax)
TukeyHSD(fit)






## Look at function data ####
rates=read.csv("16S/CH4CO2rates.csv", stringsAsFactors = F)
rates=merge(rates, meta, by="sampleID")
rates=rates[rates$lakeID.x!="CN",]
#rates=rates[!(rates$lakeID.x=="PR"&rates$treatment.x=="acidic"&rates$replicate.x==3),]
rates$treatment.x=factor(rates$treatment.x,
                         levels=c("acidic", "neutral"),
                         labels=c("acidic", "near-neutral"))
rates$lakeType=factor(rates$lakeType,
                         levels=c("acidic", "neutral", "regional"),
                         labels=c("A", "N", "R"))
p4=ggplot(rates, aes(x=lakeType, y=CH4prod, fill=treatment.x))+
  geom_boxplot()+xlab("Inoculum Source")+
  ylab(expression(paste(CH["4"], " production ", sep="")))+
  scale_fill_discrete(name="pH Environment")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.position = c(0.75, 0.8),
        legend.background = element_rect(fill="grey90"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9))
p5=ggplot(rates, aes(x=lakeType, y=CO2prod, fill=treatment.x))+
  geom_boxplot()+xlab("Inoculum Source")+
  ylab(expression(paste(CO["2"], " production ", sep="")))+
  scale_fill_discrete(name="pH Environment")+
  theme(legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.3,"cm"), 
        legend.justification = "top", 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9))
p6=plot_grid(p4, p5+guides(fill=F),
             ncol=1)
ggsave("figures/function.pdf",p6, height=4.5, width=3)

fit1=aov(CH4prod~treatment.x*lakeType, data=rates[rates$lakeID.x!="RP",])

t.test(CH4prod~lakeType, data=rates[rates$lakeType!="A",])
t.test(CO2prod~treatment.x, data=rates[rates$lakeType=="R",])

fit1sum=summary(fit1)
#Treatment effect size 
F_to_eta2(f=fit1sum[[1]]$`F value`[1], df=1, df_error = 50)
#Lake type effect size
F_to_eta2(f=fit1sum[[1]]$`F value`[2], df=2,df_error = 50)
#Interaction effect size
F_to_eta2(f=fit1sum[[1]]$`F value`[3], df=2,df_error = 50)


fit2=aov(CO2prod~treatment.x*lakeType, data=rates[rates$lakeID.x!="RP",])
fit2sum=summary(fit2)
#Treatment effect size = 0.09
F_to_eta2(f=fit2sum[[1]]$`F value`[1], df=1, df_error = 50)
#Lake type effect size = 0.23
F_to_eta2(f=fit2sum[[1]]$`F value`[2], df=2,df_error = 50)
#Interaction effect size = 0.09
F_to_eta2(f=fit2sum[[1]]$`F value`[3], df=2,df_error = 50)


