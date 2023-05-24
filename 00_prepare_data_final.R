### Prepare data
library(nadiv);library(tidyverse)
setwd("~/../Desktop/reltedness_2k22/")
geno = read.table('./DRYAD/data/genotypeTable',na.strings = "NA", header=T, sep=" ")#load data
map = read.table('./DRYAD/data/genotypeMap')#load map
pheno = read.table('./DRYAD/data/complete_pheno_data.txt',header = TRUE,sep="\t")

geno$id[geno$id=="42-f"] <- "42-f-NA"# need to rename one individual to match with pheno id
HELdata = inner_join(pheno,geno, by="id", suffix=c("",".y"))%>%
              select(-ends_with(".y"))# make total dataset

# Make the parents wild only data
parents = HELdata[grepl("-",HELdata$id),]

# Make genotype matrices for Helsinki dataset
X1 <- HELdata[,20:49512]
X1 = X1[,-which(map[,1]==12)] #remove sex chromosome from data
X1 = as.matrix(X1)
X2 = parents[,20:49512]
X2 = X2[,-which(map[,1]==12)] #remove sex chromosome from data
X2 = as.matrix(X2)