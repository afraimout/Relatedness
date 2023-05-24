#### Table S1

h2s = c(median(h2_SL_PED), median(h2_SL_SNP), median(h2_SL_SNPFS), median(h2_SL_SNPWILD),
        median(h2_BD_PED), median(h2_BD_SNP), median(h2_BD_SNPFS),
        median(h2_PL_PED), median(h2_PL_SNP), median(h2_PL_SNPFS))

Vas = c(median(VaSLPED_adj), median(VaSLSNP_adj), median(VaSLSNPFS_adj), median(VaSLSNPWILD_adj),
        median(VaBDPED_adj), median(VaBDSNP_adj), median(VaBDSNPFS_adj),
        median(VaPLPED_adj), median(VaPLSNP_adj), median(VaPLSNPFS_adj))

Vds = c(median(VdSLPED_adj), median(VdSLSNP_adj), median(VdSLSNPFS_adj), median(VdSLSNPWILD_adj),
        median(VdBDPED_adj), median(VdBDSNP_adj), median(VdBDSNPFS_adj),
        median(VdPLPED_adj), median(VdPLSNP_adj), median(VdPLSNPFS_adj))

Vps = c(median(VpSLPED_adj), median(VpSLSNP_adj), median(VpSLSNPFS_adj), median(VpSLSNPWILD_adj),
        median(VpBDPED_adj), median(VpBDSNP_adj), median(VpBDSNPFS_adj),
        median(VpPLPED_adj), median(VpPLSNP_adj), median(VpPLSNPFS_adj))

Vrs = c(median(VrSLPED_adj), median(VrSLSNP_adj), median(VrSLSNPFS_adj), median(VrSLSNPWILD_adj),
        median(VrBDPED_adj), median(VrBDSNP_adj), median(VrBDSNPFS_adj),
        median(VrPLPED_adj), median(VrPLSNP_adj), median(VrPLSNPFS_adj))

QG = data.frame(h2 = h2s,
                Va = Vas,
                Vd = Vds,
                Vp = Vps,
                Vr = Vrs)

Valist = list(VaSLPED_adj, VaSLSNP_adj, VaSLSNPFS_adj, VaSLSNPWILD_adj,
              VaBDPED_adj, VaBDSNP_adj, VaBDSNPFS_adj,
              VaPLPED_adj, VaPLSNP_adj, VaPLSNPFS_adj)
Vas.ci =lapply(Valist, function(x) HPDinterval(x))

Vdlist = list(VdSLPED_adj, VdSLSNP_adj, VdSLSNPFS_adj, VdSLSNPWILD_adj,
              VdBDPED_adj, VdBDSNP_adj, VdBDSNPFS_adj,
              VdPLPED_adj, VdPLSNP_adj, VdPLSNPFS_adj)
Vds.ci =lapply(Vdlist, function(x) HPDinterval(x))

Vrlist = list(VrSLPED_adj, VrSLSNP_adj, VrSLSNPFS_adj, VrSLSNPWILD_adj,
              VrBDPED_adj, VrBDSNP_adj, VrBDSNPFS_adj,
              VrPLPED_adj, VrPLSNP_adj, VrPLSNPFS_adj)
Vrs.ci =lapply(Vrlist, function(x) HPDinterval(x))

Vplist = list(VpSLPED_adj, VpSLSNP_adj, VpSLSNPFS_adj, VpSLSNPWILD_adj,
              VpBDPED_adj, VpBDSNP_adj, VpBDSNPFS_adj,
              VpPLPED_adj, VpPLSNP_adj, VpPLSNPFS_adj)
Vps.ci =lapply(Vplist, function(x) HPDinterval(x))

h2list = list(h2_SL_PED, h2_SL_SNP, h2_SL_SNPFS, h2_SL_SNPWILD,
              h2_BD_PED, h2_BD_SNP, h2_BD_SNPFS,
              h2_PL_PED, h2_PL_SNP, h2_PL_SNPFS)
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

QGtable$type = factor(c("SL_PED", "SL_HS","SL_FS","SL_Wild",
                        "BD_PED", "BD_HS","BD_FS",
                        "PL_PED", "PL_HS","PL_FS"))

#### gather
QG1 =  QGtable %>% gather(varcomp, median, c("Va","Vd","Vr","Vp","h2"))
QG2 =  QGtable %>% gather(varcomp, ci.low, c("Va.LCI","Vd.LCI","Vr.LCI","Vp.LCI","h2.LCI"))
QG3 =  QGtable %>% gather(varcomp, ci.up, c("Va.UCI","Vd.UCI","Vr.UCI","Vp.UCI","h2.UCI"))

QG = data.frame(trait = substr(QG1$type, 1,2),
                cross = substr(QG1$type, 4,7),
                varcomp = factor(QG1$varcomp),
                median = round(QG1$median,3),
                ci.low = as.numeric(formatC(QG2$ci.low, format = "e", digits = 2)),
                ci.up = round(QG3$ci.up,3))


QG$type = factor(paste(QG$trait,QG$cross,sep="_"),levels = c("SL_PED","SL_HS","SL_FS","SL_Wild",
                                                             "BD_PED","BD_HS","BD_FS", 
                                                             "PL_PED","PL_HS","PL_FS"))


Table_S1a=QG[QG$cross=="PED",]
knitr::kable(Table_S1a)  
Table_S1b=QG[QG$cross=="HS",]
knitr::kable(Table_S1b)  
