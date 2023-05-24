### Calculate the Dk statistic from Legarra (2015) to correct PRM variance estimates
setwd("~/../Desktop/reltedness_2k22/DRYAD/")
load("./data/DataAndGRMsHelsinki.RData")
library(MCMCglmm)
library(nadiv)
path = paste0(getwd(),"/results/Animal_models/")
filenames = list.files(path, pattern = ".Rdata")

for (i in 1:length(filenames)){
  load(paste0(path,filenames[i]))
}

# load the PRM animal models
load("./results/Animal_models/model_SL_helsinki_PED.Rdata")
load("./results/Animal_models/model_BD_helsinki_PED.Rdata")
load("./results/Animal_models/model_PL_helsinki_PED.Rdata")

# load the GRM animal models
load("./results/Animal_models/model_SL_helsinki_SNP.Rdata")
load("./results/Animal_models/model_BD_helsinki_SNP.Rdata")
load("./results/Animal_models/model_PL_helsinki_SNP.Rdata")

SLPED = lapply(mod_SL_ped, function(m) m$VCV)
BDPED = lapply(mod_BD_ped, function(m) m$VCV)
PLPED = lapply(mod_PL_ped, function(m) m$VCV)

SLSNP = model_SL$VCV
BDSNP = model_BD$VCV
PLSNP = model_PL$VCV

SLSNPFS = model_SL_FS$VCV
BDSNPFS = model_BD_FS$VCV
PLSNPFS = model_PL_FS$VCV

SLSNPWILD = model_SL_helsinki_parents_GRM$VCV

### QG estimates
#SL
VaSLPED = SLPED[[1]][, "animal"]
VaSLSNP = SLSNP[, "id"]
VaSLSNPFS = SLSNPFS[, "id"]
VaSLSNPWILD = SLSNPWILD[, "id"]

VdSLPED = SLPED[[1]][, "dom"]
VdSLSNP = SLSNP[, "dom"]
VdSLSNPFS = SLSNPFS[, "dom"]
VdSLSNPWILD = SLSNPWILD[, "dom"]

VrSLPED = SLPED[[1]][, "units"]
VrSLSNP = SLSNP[, "units"]
VrSLSNPFS = SLSNPFS[, "units"]
VrSLSNPWILD = SLSNPWILD[, "units"]

#BD
VaBDPED = BDPED[[1]][, "animal"]
VaBDSNP = BDSNP[, "id"]
VaBDSNPFS = BDSNPFS[, "id"]

VdBDPED = BDPED[[1]][, "dom"]
VdBDSNP = BDSNP[, "dom"]
VdBDSNPFS = BDSNPFS[, "dom"]

VrBDPED = BDPED[[1]][, "units"]
VrBDSNP = BDSNP[, "units"]
VrBDSNPFS = BDSNPFS[, "units"]

#PL
VaPLPED = PLPED[[1]][, "animal"]
VaPLSNP = PLSNP[, "id"]
VaPLSNPFS = PLSNPFS[, "id"]

VdPLPED = PLPED[[1]][, "dom"]
VdPLSNP = PLSNP[, "dom"]
VdPLSNPFS = PLSNPFS[, "dom"]

VrPLPED = PLPED[[1]][, "units"]
VrPLSNP = PLSNP[, "units"]
VrPLSNPFS = PLSNPFS[, "units"]

# Pedigree relationship matrices
pedHEL = HELdata[,4:6] #pedigree for Helsinki
pedHEL[pedHEL == 0] <- NA #missing parents must be coded as NAs for nadiv
HELdata$dom = HELdata$animal #adding extra id column as dom to estimate variance dominance

HelD = makeD(pedHEL, invertD = FALSE)
HelA = makeA(pedHEL)

HelA = as.matrix(HelA)
diag(HelA)= 1
HelD = as.matrix(HelD$D)
HelR = model_SL$ZR

#Calculate Dk
n = dim(pedHEL)[1]

DK_A = ((sum(diag(HelA))/n))-mean(HelA)
DK_D = ((sum(diag(HelD))/n))-mean(HelD)
DK_R = ((sum(diag(HelR))/n))-mean(HelR)
DK_GA = ((sum(diag(GA_Hel))/n))-mean(GA_Hel)
DK_GD = ((sum(diag(GD_Hel))/n))-mean(GD_Hel)

#Apply correction
VaSLPED_adj = VaSLPED * DK_A
VdSLPED_adj = VdSLPED * DK_D
VrSLPED_adj = VrSLPED * DK_R

VaSLSNP_adj = VaSLSNP * DK_GA
VdSLSNP_adj = VdSLSNP * DK_GD
VrSLSNP_adj = VrSLSNP * DK_R

VaSLSNPFS_adj = VaSLSNPFS * DK_GA
VdSLSNPFS_adj = VdSLSNPFS * DK_GD
VrSLSNPFS_adj = VrSLSNPFS * DK_R

VaSLSNPWILD_adj = VaSLSNPWILD * DK_GA
VdSLSNPWILD_adj = VdSLSNPWILD * DK_GD
VrSLSNPWILD_adj = VrSLSNPWILD * DK_R

VpSLPED_adj = VaSLPED_adj + VdSLPED_adj + VrSLPED_adj
VpSLSNP_adj = VaSLSNP_adj + VdSLSNP_adj + VrSLSNP_adj
VpSLSNPFS_adj = VaSLSNPFS_adj + VdSLSNPFS_adj + VrSLSNPFS_adj
VpSLSNPWILD_adj = VaSLSNPWILD_adj + VdSLSNPWILD_adj + VrSLSNPWILD_adj

h2_SL_PED = VaSLPED_adj / (VpSLPED_adj)
h2_SL_SNP = VaSLSNP_adj / (VpSLSNP_adj)
h2_SL_SNPFS = VaSLSNPFS_adj / (VpSLSNPFS_adj)
h2_SL_SNPWILD = VaSLSNPWILD_adj / (VpSLSNPWILD_adj)

#BD
VaBDPED_adj = VaBDPED * DK_A
VdBDPED_adj = VdBDPED * DK_D
VrBDPED_adj = VrBDPED * DK_R

VaBDSNP_adj = VaBDSNP * DK_GA
VdBDSNP_adj = VdBDSNP * DK_GD
VrBDSNP_adj = VrBDSNP * DK_R

VaBDSNPFS_adj = VaBDSNPFS * DK_GA
VdBDSNPFS_adj = VdBDSNPFS * DK_GD
VrBDSNPFS_adj = VrBDSNPFS * DK_R

VpBDPED_adj = VaBDPED_adj + VdBDPED_adj + VrBDPED_adj
VpBDSNP_adj = VaBDSNP_adj + VdBDSNP_adj + VrBDSNP_adj
VpBDSNPFS_adj = VaBDSNPFS_adj + VdBDSNPFS_adj + VrBDSNPFS_adj

h2_BD_PED = VaBDPED_adj / (VaBDPED_adj + VdBDPED_adj + VrBDPED_adj)
h2_BD_SNP = VaBDSNP_adj / (VaBDSNP_adj + VdBDSNP_adj + VrBDSNP_adj)
h2_BD_SNPFS = VaBDSNPFS_adj / (VaBDSNPFS_adj + VdBDSNPFS_adj + VrBDSNPFS_adj)

#PL
VaPLPED_adj = VaPLPED * DK_A
VdPLPED_adj = VdPLPED * DK_D
VrPLPED_adj = VrPLPED * DK_R

VaPLSNP_adj = VaPLSNP * DK_GA
VdPLSNP_adj = VdPLSNP * DK_GD
VrPLSNP_adj = VrPLSNP * DK_R

VaPLSNPFS_adj = VaPLSNPFS * DK_GA
VdPLSNPFS_adj = VdPLSNPFS * DK_GD
VrPLSNPFS_adj = VrPLSNPFS * DK_R

VpPLPED_adj = VaPLPED_adj + VdPLPED_adj + VrPLPED_adj
VpPLSNP_adj = VaPLSNP_adj + VdPLSNP_adj + VrPLSNP_adj
VpPLSNPFS_adj = VaPLSNPFS_adj + VdPLSNPFS_adj + VrPLSNPFS_adj

h2_PL_PED = VaPLPED_adj / (VaPLPED_adj + VdPLPED_adj + VrPLPED_adj)
h2_PL_SNP = VaPLSNP_adj / (VaPLSNP_adj + VdPLSNP_adj + VrPLSNP_adj)
h2_PL_SNPFS = VaPLSNPFS_adj / (VaPLSNPFS_adj + VdPLSNPFS_adj + VrPLSNPFS_adj)

h2s = c(median(h2_SL_SNP), median(h2_SL_SNPFS), median(h2_SL_SNPWILD),
        median(h2_BD_SNP), median(h2_BD_SNPFS),
        median(h2_PL_SNP), median(h2_PL_SNPFS))

Vas = c(median(VaSLSNP_adj), median(VaSLSNPFS_adj), median(VaSLSNPWILD_adj),
        median(VaBDSNP_adj), median(VaBDSNPFS_adj),
        median(VaPLSNP_adj), median(VaPLSNPFS_adj))

Vds = c(median(VdSLSNP_adj), median(VdSLSNPFS_adj), median(VdSLSNPWILD_adj),
        median(VdBDSNP_adj), median(VdBDSNPFS_adj),
        median(VdPLSNP_adj), median(VdPLSNPFS_adj))

Vps = c(median(VpSLSNP_adj), median(VpSLSNPFS_adj), median(VpSLSNPWILD_adj),
        median(VpBDSNP_adj), median(VpBDSNPFS_adj),
        median(VpPLSNP_adj), median(VpPLSNPFS_adj))

Vrs = c(median(VrSLSNP_adj), median(VrSLSNPFS_adj), median(VrSLSNPWILD_adj),
        median(VrBDSNP_adj), median(VrBDSNPFS_adj),
        median(VrPLSNP_adj), median(VrPLSNPFS_adj))

QG = data.frame(h2 = h2s,
                Va = Vas,
                Vd = Vds,
                Vp = Vps,
                Vr = Vrs)

Valist = list(VaSLSNP_adj, VaSLSNPFS_adj, VaSLSNPWILD_adj,
              VaBDSNP_adj, VaBDSNPFS_adj,
              VaPLSNP_adj, VaPLSNPFS_adj)
Vas.ci =lapply(Valist, function(x) HPDinterval(x))

Vdlist = list(VdSLSNP_adj, VdSLSNPFS_adj, VdSLSNPWILD_adj,
              VdBDSNP_adj, VdBDSNPFS_adj,
              VdPLSNP_adj, VdPLSNPFS_adj)
Vds.ci =lapply(Vdlist, function(x) HPDinterval(x))

Vrlist = list(VrSLSNP_adj, VrSLSNPFS_adj, VrSLSNPWILD_adj,
              VrBDSNP_adj, VrBDSNPFS_adj,
              VrPLSNP_adj, VrPLSNPFS_adj)
Vrs.ci =lapply(Vrlist, function(x) HPDinterval(x))

Vplist = list(VpSLSNP_adj, VpSLSNPFS_adj, VpSLSNPWILD_adj,
              VpBDSNP_adj, VpBDSNPFS_adj,
              VpPLSNP_adj, VpPLSNPFS_adj)
Vps.ci =lapply(Vplist, function(x) HPDinterval(x))

h2list = list(h2_SL_SNP, h2_SL_SNPFS, h2_SL_SNPWILD,
              h2_BD_SNP, h2_BD_SNPFS,
              h2_PL_SNP, h2_PL_SNPFS)
h2s.ci =lapply(h2list, function(x) HPDinterval(x))

QGtable = data.frame(Va = NA,
                     Va.LCI = NA,
                     Va.UCI = NA,
                     Vd = NA,
                     Vd.LCI = NA,
                     Vd.UCI = NA,
                     Vr = NA,
                     Vr.LCI = NA,
                     Vr.UCI = NA,
                     Vp = NA,
                     Vp.LCI = NA,
                     Vp.UCI = NA,
                     h2 = NA,
                     h2.LCI = NA,
                     h2.UCI = NA)

for (i in 1:length(Vas)) {
  QGtable[i,1] <- unlist(Vas[[i]])
  QGtable[i,2] <- unlist(Vas.ci[[i]])[,1]
  QGtable[i,3] <- unlist(Vas.ci[[i]])[,2]
  QGtable[i,4] <- unlist(Vds[[i]])
  QGtable[i,5] <- unlist(Vds.ci[[i]])[,1]
  QGtable[i,6] <- unlist(Vds.ci[[i]])[,2]
  QGtable[i,7] <- unlist(Vrs[[i]])
  QGtable[i,8] <- unlist(Vrs.ci[[i]])[,1]
  QGtable[i,9] <- unlist(Vrs.ci[[i]])[,2]
  QGtable[i,10] <- unlist(Vps[[i]])
  QGtable[i,11] <- unlist(Vps.ci[[i]])[,1]
  QGtable[i,12] <- unlist(Vps.ci[[i]])[,2]
  QGtable[i,13] <- unlist(h2s[[i]])
  QGtable[i,14] <- unlist(h2s.ci[[i]])[,1]
  QGtable[i,15] <- unlist(h2s.ci[[i]])[,2]
}

####### plot QG params
library(tidyverse)
library(PNWColors)
pal = pnw_palette(name="Sunset",n=7,type="discrete")

QGtable$type = factor(c("SL_HS","SL_FS","SL_Wild","BD_HS","BD_FS","PL_HS","PL_FS"))

#### gather
QG1 =  QGtable %>% gather(varcomp, median, c("Va","Vd","Vr","Vp","h2"))
QG2 =  QGtable %>% gather(varcomp, ci.low, c("Va.LCI","Vd.LCI","Vr.LCI","Vp.LCI","h2.LCI"))
QG3 =  QGtable %>% gather(varcomp, ci.up, c("Va.UCI","Vd.UCI","Vr.UCI","Vp.UCI","h2.UCI"))

QG = data.frame(trait = substr(QG1$type, 1,2),
                cross = substr(QG1$type, 4,7),
                varcomp = factor(QG1$varcomp),
                median = QG1$median,
                ci.low = QG2$ci.low,
                ci.up = QG3$ci.up)

QG$type = factor(paste(QG$trait,QG$cross,sep="_"),levels = c("SL_HS","SL_FS","SL_Wild",
                                                               "BD_HS","BD_FS", 
                                                               "PL_HS","PL_FS"))
Fig1=
  QG %>%
  ggplot(aes(x = type, y = median, color = varcomp,shape=cross))+
  geom_point(size=3) +
  geom_errorbar(aes(ymin = ci.low, ymax = ci.up,width=.02), linewidth=1) +
  theme_bw() +
  facet_wrap(trait~varcomp,ncol = 5,scales="free")+
  guides(color="none")+
  scale_shape_manual(values = c(15, 16, 17))+
  scale_color_manual(values=pal[1:5],
                     aesthetics=c("fill","color")) +
  theme_bw() +
  ylab(expression(paste("")))+
  xlab(expression(paste("")))+
  theme(aspect.ratio=1,
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank())+
  theme(axis.text.x=element_text(angle=45, hjust=1,face = "bold",size=6))
