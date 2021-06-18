# Functions for indicator analyses


#### Get indicator ZOTUs associated with pH environments
#data=indvalM.pH$sign
extractIndicatorZOTUs.pH=function(data){
  data=data[!is.na(data$p.value),]
  aZOTUs=data.frame(pH="acidic", ZOTUs=rownames(data[data$s.1==1 & data$p.value<=0.05,]))
  nZOTUs=data.frame(pH="near-neutral", ZOTUs=rownames(data[data$s.2==1 & data$p.value<=0.05,]))
  return(rbind(aZOTUs, nZOTUs))
}

#### Get indicator ZOTUs associated with inoculum source
#data=indvalM.inoc$sign
extractIndicatorZOTUs.inoc=function(data){
  data=data[!is.na(data$p.value),]
  aZOTUs=data.frame(source="acidic", ZOTUs=rownames(data[data$s.1==1 & data$p.value<=0.05,]))
  nZOTUs=data.frame(source="near-neutral", ZOTUs=rownames(data[data$s.2==1 & data$p.value<=0.05,]))
  return(rbind(aZOTUs, nZOTUs))
}

#### Generate heat map for MCC
genHeatmap.MCC=function(ind, label){ #ind=indvalM.pH.out
  colnames(ind)=c("treatment", "V1")
  
  ind=merge(ind, tax, by="V1")

  ind2=ind %>% separate(col="V2", into=c("d", "p", "c", "o", "f", "g"), sep = ",[dpcofg]:")
  ind2$sum=1
  ind2$p=gsub("\\([01][.][0-9]{4}\\)", "", ind2$p)
  ind2$o=gsub("\\([01][.][0-9]{4}\\)", "", ind2$o)
  ind2$g=gsub("\\([01][.][0-9]{4}\\)", "", ind2$g)
  ind3=aggregate(sum~g+treatment, data=ind2, FUN="sum")
  
  p1=ggplot(ind3, aes(x=treatment, y=g))+
    geom_tile(aes(fill=sum))+
    geom_text(aes(label=sum), color="white")+
    scale_fill_gradient(low="grey75", high="black")+
    guides(fill=F)+
    xlab(label)+ylab("Genus")
  return(p1)
}


#### Generate heat map for MCC
genHeatmap.nonMCC=function(ind, label){ #ind=indvalM.pH.out
  colnames(ind)=c("treatment", "V1")
  
  ind=merge(ind, tax, by="V1")
  
  ind2=ind %>% separate(col="V2", into=c("d", "p", "c", "o", "f", "g"), sep = ",[dpcofg]:")
  ind2$sum=1
  ind2$p=gsub("\\([01][.][0-9]{4}\\)", "", ind2$p)
  ind2$o=gsub("\\([01][.][0-9]{4}\\)", "", ind2$o)
  ind2$g=gsub("\\([01][.][0-9]{4}\\)", "", ind2$g)
  ind3=aggregate(sum~p+treatment, data=ind2, FUN="sum")
  
  p1=ggplot(ind3, aes(x=treatment, y=p))+
    geom_tile(aes(fill=sum))+
    geom_text(aes(label=sum), color="white")+
    scale_fill_gradient(low="grey75", high="black")+
    guides(fill=F)+
    xlab(label)+ylab("Phylum")
  return(p1)
}
