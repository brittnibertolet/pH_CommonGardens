#### tree analysis

library(ape)
library(dplyr)
library(rlang)
source("code/delta_statistic_functions.R")


#### Look at phylogenetic signal of methanogen across inoculum ####
#Read in tree and here we take 1% of the 1% quantile to fill in the null branches
tree=read.tree("trees/indicator.M.inoc.tree")
tree=multi2di(tree)
tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1

#Read in trait information
indM.inoc=read.csv("trees/indicator.M.inoc2.csv", stringsAsFactors = F)

#Get trait information 
orderTree=data.frame(ZOTUs=tree$tip.label)
orderTree=left_join(orderTree,indM.inoc, by="ZOTUs")
trait=orderTree$source

#Calculate delta 
deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100) #2.63

#Get pvalue by comparing against deltas generated from random trait data
random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>deltaA)/length(random_delta)
boxplot(random_delta)
abline(h=deltaA,col="red")

#### Look at phylogenetic signal of non-methanogen across inoculum ####
#Read in tree and here we take 1% of the 1% quantile to fill in the null branches
tree=read.tree("trees/indicator.B.inoc.tree")
tree=multi2di(tree)
tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1
tree=indM.inoc.tree

#Read in trait information
indM.inoc=read.csv("trees/indicator.B.inoc2.csv", stringsAsFactors = F)
indM.inoc=indM.inoc[!duplicated(indM.inoc$ZOTUs),]

#Get trait information 
orderTree=data.frame(ZOTUs=tree$tip.label)
orderTree=left_join(orderTree,indM.inoc, by="ZOTUs")
trait=orderTree$source

#Calculate delta
deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100) #2.22

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>deltaA)/length(random_delta)
boxplot(random_delta)
abline(h=deltaA,col="red")

#### Look at phylogenetic signal of methanogen across pH #### 
#Read in tree and here we take 1% of the 1% quantile to fill in the null branches
tree=read.tree("trees/indicator.M.pH.tree")
tree=multi2di(tree)
tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1

#Read in trait information
indM.inoc=read.csv("trees/indicator.M.pH2.csv", stringsAsFactors = F)
indM.inoc=indM.inoc[!duplicated(indM.inoc$ZOTUs),]

#Get trait information 
orderTree=data.frame(ZOTUs=tree$tip.label)
orderTree=left_join(orderTree,indM.inoc, by="ZOTUs")
trait=orderTree$pH

#Calculate delta
deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100) #3.96

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>deltaA)/length(random_delta) # p = 0
boxplot(random_delta, ylim=c(0,4.0))
abline(h=deltaA,col="red")


#### Look at phylogenetic signal of nonmethanogen across pH ####
#Read in tree and here we take 1% of the 1% quantile to fill in the null branches
tree=read.tree("trees/indicator.B.pH.tree")
tree=multi2di(tree)
tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1
tree=indM.inoc.tree

#Read in trait information
indM.inoc=read.csv("trees/indicator.B.pH2.csv", stringsAsFactors = F)
indM.inoc=indM.inoc[!duplicated(indM.inoc$ZOTUs),]

#Get trait information 
orderTree=data.frame(ZOTUs=tree$tip.label)
orderTree=left_join(orderTree,indM.inoc, by="ZOTUs")
trait=orderTree$pH

#Calculate delta
deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100) #3.86

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>deltaA)/length(random_delta) # 0
boxplot(random_delta, ylim=c(0,4))
abline(h=deltaA,col="red")
