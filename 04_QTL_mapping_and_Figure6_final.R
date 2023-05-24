######### singlemapping 4 way for 3 crosses ##########
# Script to perform singlamapping 4-way for the three 9sp crosses
setwd("~/../Desktop/reltedness_2k22/DRYAD/results/QTL_mapping/")
########## prepare data ##############
load("data_3crosses.RData")
source('singlemapping_fourway.R')
source("plotting_functions.R")
library(data.table);library(ggplot2)
library(foreach);library(parallel)

# id data for genotypes
IND <- data.table(IND)
setnames(IND, colnames(IND), c("cross", "ID", "sex"))

# get covariates (sex and SL)
sexR <- sex[which(IND[,1]=='HR')]
sexP <- sex[which(IND[,1]=='HP')]
sexB <- sex[which(IND[,1]=='HB')]

# remove sex linked markers from complexity reduced data
IDsexR <- which(abs(cor(X_cl[[1]],sexR))>0.95)
IDsexP <- which(abs(cor(X_cl[[4]],sexP))>0.95)
IDsexB <- which(abs(cor(X_cl[[7]],sexB))>0.95)

XR <- cbind(X_cl[[1]][,-IDsexR],X_cl[[2]][,-IDsexR],X_cl[[3]][,-IDsexR])
XP <- cbind(X_cl[[4]][,-IDsexR],X_cl[[5]][,-IDsexR],X_cl[[6]][,-IDsexR])
XB <- cbind(X_cl[[7]][,-IDsexR],X_cl[[8]][,-IDsexR],X_cl[[9]][,-IDsexR])

# get map for complexity reduced data
cluster_summary <- data.table(CL_data_9sp_codom[[1]])
map_cl <- as.matrix(cluster_summary[,.(Chr, Pos)][-IDsexR,])

# get data for each cross, 1 is male, 2 is female, 3 is dominance effect
p <- dim(XR)[2]/3
XR_1 <- XR[,1:p]
XR_2 <- XR[,(p+1):(2*p)]
XR_3 <- XR[,(2*p+1):(3*p)]

XP_1 <- XP[,1:p]
XP_2 <- XP[,(p+1):(2*p)]
XP_3 <- XP[,(2*p+1):(3*p)]

XB_1 <- XB[,1:p]
XB_2 <- XB[,(p+1):(2*p)]
XB_3 <- XB[,(2*p+1):(3*p)]

########## prepare phenotypes ###############
load("QTL_relatedness_pheno_data.RData")

phe_R = as.data.table(HELRYT)
phe_P = as.data.table(HELPYO)
phe_B = as.data.table(HELBYN)

phe_R[,ID := gsub("-", "-", phe_R[,ID])] 
phe_P[,ID := gsub("-", "-", phe_P[,ID])] 
phe_B[,ID := gsub("–", "-", phe_B[,ID])] 

## match with genotype data, IND data needs to be prepared first, see Start_here
phe_R <- phe_R[match(as.character(IND[cross=="HR",ID]),phe_R[,ID]),]
phe_B <- phe_B[match(as.character(IND[cross=="HB",ID]),as.character(phe_B[,ID])),]
phe_P <- phe_P[match(as.character(IND[cross=="HP",ID]),phe_P[,ID]),]

########## singlemapping 4-way ##############
# Singelmapping 4-way, "phe_" is a matrix with individuals in rows and phenotypes as columns, matched with IND[,2]. These are sorted out in "prepare_phenotypes.R"
cores <- 1
nrep <- 10000 # number of repetitions for permutation, 0 means no permutation and you only get uncorrected pvalues

### run
single4wayR <- apply(phe_R[,3:5], 2, function(y){
  list(singlemapping_fourway(X=XR_1[!is.na(y),],X1=XR_2[!is.na(y),],X2=XR_3[!is.na(y),],y=y[!is.na(y)],covariate=sexR[!is.na(y)],num_perm=nrep,core=cores),
       singlemapping_fourway(X=XR_2[!is.na(y),],X1=XR_1[!is.na(y),],X2=XR_3[!is.na(y),],y=y[!is.na(y)],covariate=sexR[!is.na(y)],num_perm=nrep,core=cores),
       singlemapping_fourway(X=XR_3[!is.na(y),],X1=XR_1[!is.na(y),],X2=XR_2[!is.na(y),],y=y[!is.na(y)],covariate=sexR[!is.na(y)],num_perm=nrep,core=cores))  
  
})
#save(single4wayR, file = "single4wayR_SLBDPL_final.Rdata")

single4wayB <- apply(phe_B[,3:5], 2, function(y){
  list(singlemapping_fourway(X=XB_1[!is.na(y),],X1=XB_2[!is.na(y),],X2=XB_3[!is.na(y),],y=y[!is.na(y)],covariate=sexB[!is.na(y)],num_perm=nrep,core=cores),
       singlemapping_fourway(X=XB_2[!is.na(y),],X1=XB_1[!is.na(y),],X2=XB_3[!is.na(y),],y=y[!is.na(y)],covariate=sexB[!is.na(y)],num_perm=nrep,core=cores),
       singlemapping_fourway(X=XB_3[!is.na(y),],X1=XB_1[!is.na(y),],X2=XB_2[!is.na(y),],y=y[!is.na(y)],covariate=sexB[!is.na(y)],num_perm=nrep,core=cores))  
  
})
#save(single4wayB, file = "single4wayB_SLBDPL_final.Rdata")

single4wayP <- apply(phe_P[,3:5], 2, function(y){
  list(singlemapping_fourway(X=XP_1[!is.na(y),],X1=XP_2[!is.na(y),],X2=XP_3[!is.na(y),],y=y[!is.na(y)],covariate=sexP[!is.na(y)],num_perm=nrep,core=cores),
       singlemapping_fourway(X=XP_2[!is.na(y),],X1=XP_1[!is.na(y),],X2=XP_3[!is.na(y),],y=y[!is.na(y)],covariate=sexP[!is.na(y)],num_perm=nrep,core=cores),
       singlemapping_fourway(X=XP_3[!is.na(y),],X1=XP_1[!is.na(y),],X2=XP_2[!is.na(y),],y=y[!is.na(y)],covariate=sexP[!is.na(y)],num_perm=nrep,core=cores))  
  
})
#save(single4wayP, file = "single4wayP_SLBDPL_final.Rdata")

########## manhattan plots for all three crosses ###########
# Load this if you do not want to run the above code
load("single4wayR_SLBDPL_final.Rdata")
load("single4wayB_SLBDPL_final.Rdata")
load("single4wayP_SLBDPL_final.Rdata")

## some data needed for the plots
MAP2 <- data.frame(map_cl)
MAP2$Chr <- paste("chr",as.character(MAP2$Chr),sep="")

genome=data.frame(paste("chr",as.character(as.character(1:21)),sep=""), tapply(MAP2[,2], as.numeric(gsub("chr", "",  MAP2[,1])), max))

colnames(genome) <- c("V1","V2")

## some useful colors for chromosomes
col_vector <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", 
                "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", 
                "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
                "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", 
                "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", 
                "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", 
                "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", 
                "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
                "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
                "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
                "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", 
                "#CCEBC5", "#FFED6F")

## use this to find the order of colors that gives nice separation between adjacent pairs of chromosomes
col_manhattan <- sample(col_vector)

################## 
## get P values ##
################## 

##For PCs
phenotypes <- phe_R
trait_names <- colnames(phenotypes)[3:5] ## early late growth

HR_p_val <- lapply(1:length(single4wayR), function(i){
  lapply(1:3, function(x){
    temp <- single4wayR[[i]][[x]]$pval.corr ## use this and comment out the below to use corrected p.values (when nrep =>1000)
    #temp <- single4wayR_bin[[i]][[x]]$pval
    temp[temp == 0] <- 1/1000
    temp
  })
  
})


HB_p_val <- lapply(1:length(single4wayB), function(i){
  lapply(1:3, function(x){
    temp <- single4wayB[[i]][[x]]$pval.corr 
    #temp <- single4wayB_bin[[i]][[x]]$pval
    temp[temp == 0] <- 1/1000
    temp
  })
  
})

names(HB_p_val) <- trait_names

HP_p_val <- lapply(1:length(single4wayP), function(i){
  lapply(1:3, function(x){
    temp <- single4wayP[[i]][[x]]$pval.corr
    #temp <- single4wayP_bin[[i]][[x]]$pval
    temp[temp == 0] <- 1/1000
    temp
  })
  
})

names(HP_p_val) <- trait_names

#### get necessary data for manhattan plots 

#HELBYN
dt_HB <- lapply(1:length(HB_p_val), function(i){
  lapply(1:3, function(x){
    pvalues=HB_p_val[[i]][[x]]
    plotManhattan_ggplot(bedfile=MAP2,pvalues=pvalues,genome=genome,col=col_manhattan,cex=0.75, space = 0.01, lwd=2, yrange = c(-0.1, 3.2), type="h")
  })
})


dt_HB <- do.call(rbind, unlist(lapply(dt_HB, function(x) lapply(x, function(y) y[[1]])), recursive = FALSE))

Labels <- expand.grid(data = c("M","F","D"), trait = trait_names)
dt_HB[,cross := rep("HB", each=241)]
dt_HB[,trait := rep(Labels$trait, each=241)]
dt_HB[,data := rep(Labels$data, each=241)]

#HELPYO
dt_HP <- lapply(1:length(HP_p_val), function(i){
  lapply(1:3, function(x){
    pvalues=HP_p_val[[i]][[x]]
    plotManhattan_ggplot(bedfile=MAP2,pvalues=pvalues,genome=genome,col=col_manhattan,cex=0.75, space = 0.01, lwd=2, yrange = c(-0.1, 3.2), type="h")
  })
})

dt_HP <- do.call(rbind, unlist(lapply(dt_HP, function(x) lapply(x, function(y) y[[1]])), recursive = FALSE))

Labels <- expand.grid(data = c("M","F","D"), trait = trait_names)

dt_HP[,cross := rep("HP", each=241)]
dt_HP[,trait := rep(Labels$trait, each=241)]
dt_HP[,data := rep(Labels$data, each=241)]

#HELRYT
dt_HR <- lapply(1:length(HR_p_val), function(i){
  lapply(1:3, function(x){
    pvalues=HR_p_val[[i]][[x]]
    plotManhattan_ggplot(bedfile=MAP2,pvalues=pvalues,genome=genome,col=col_manhattan,cex=0.75, space = 0.01, lwd=2, yrange = c(-0.1, 3.2), type="h")
  })
})

dt_HR <- do.call(rbind, unlist(lapply(dt_HR, function(x) lapply(x, function(y) y[[1]])), recursive = FALSE))

Labels <- expand.grid(data = c("M","F","D"), trait = trait_names)

dt_HR[,cross := rep("HR", each=241)]
dt_HR[,trait := rep(Labels$trait, each=241)]
dt_HR[,data := rep(Labels$data, each=241)]

### hera are the plots
space = 0.01
chromoffsets = chromOffsets(genome, space)
chromcenters = (chromoffsets[, 3] + chromoffsets[, 4])/2

MP_HB <- ggplot(dt_HB, aes(V1, V2, col=col))+
  #geom_segment() +
  ggalt::geom_lollipop(point.size=0.1) +
  facet_grid(data~trait) +
  scale_x_continuous("Chromosome", breaks=chromcenters, labels=labels) +
  theme_bw() +
  theme(#aspect.ratio=0.25,
    axis.text.x = element_text(colour = c("black",  "#00000000")),
    legend.position = "none",
    # strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    #axis.ticks.x = element_line(),
    #axis.line.x = element_line(colour = "black", size = .8, linetype = "solid"),
    #legend.title = element_blank(),
    #legend.background = element_rect(fill="white",color=NA),
    panel.background = element_rect(fill = "white", colour =NA, size = 0.5, linetype = "solid")#,
    #plot.background = element_rect(fill = "#F3F1E7",color=NA)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  ylab("-log10(P)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold",size = 5))



MP_HP <- ggplot(dt_HP, aes(V1, V2, col=col))+
  #geom_segment() +
  ggalt::geom_lollipop(point.size=0.1) +
  facet_grid(data~trait) +
  scale_x_continuous("Chromosome", breaks=chromcenters, labels=labels) +
  theme_bw() +
  theme(#aspect.ratio=0.25,
    axis.text.x = element_text(colour = c("black",  "#00000000")),
    legend.position = "none",
    # strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    #axis.ticks.x = element_line(),
    #axis.line.x = element_line(colour = "black", size = .8, linetype = "solid"),
    #legend.title = element_blank(),
    #legend.background = element_rect(fill="white",color=NA),
    panel.background = element_rect(fill = "white", colour =NA, size = 0.5, linetype = "solid")#,
    #plot.background = element_rect(fill = "#F3F1E7",color=NA)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  ylab("-log10(P)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold",size = 5))


MP_HR <- ggplot(dt_HR, aes(V1, V2, col=col))+
  #geom_segment() +
  ggalt::geom_lollipop(point.size=0.1) +
  facet_grid(data~trait) +
  scale_x_continuous("Chromosome", breaks=chromcenters, labels=labels) +
  theme_bw() +
  theme(#aspect.ratio=0.25,
    axis.text.x = element_text(colour = c("black",  "#00000000")),
    legend.position = "none",
    # strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    #axis.ticks.x = element_line(),
    #axis.line.x = element_line(colour = "black", size = .8, linetype = "solid"),
    #legend.title = element_blank(),
    #legend.background = element_rect(fill="white",color=NA),
    panel.background = element_rect(fill = "white", colour =NA, size = 0.5, linetype = "solid")#,
    #plot.background = element_rect(fill = "#F3F1E7",color=NA)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  ylab("-log10(P)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold",size = 5))


multiplot(MP_HB+ggtitle("A | HEL x BYN")+ylim(0,3), MP_HR+ggtitle("C | HEL x RYT")+ylim(0,3), MP_HP+ggtitle("C | HEL x PYO")+ylim(0,3))

