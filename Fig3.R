library(MCMCglmm)
setwd("./DRYAD/results/Simulations/results/migration/")
path=getwd()
filenames = list.files(path, pattern=".RData")

for (i in 1:length(filenames)){
  load(filenames[i])
}

mod_simu_R1 = mod_simu_R01WILD500_P05_M0008
mod_simu_R2 = mod_simu_R02WILD500_P05_M0008
mod_simu_R3 = mod_simu_R03WILD500_P05_M0008
mod_simu_R4 = mod_simu_R04WILD500_P05_M0008
mod_simu_R5 = mod_simu_R05WILD500_P05_M0008
mod_simu_R6 = mod_simu_R06WILD500_P05_M0008
mod_simu_R7 = mod_simu_R07WILD500_P05_M0008
mod_simu_R8 = mod_simu_R08WILD500_P05_M0008
mod_simu_R9 = mod_simu_R09WILD500_P05_M0008
mod_simu_R10 = mod_simu_R10WILD500_P05_M0008

mod_simu_R11 = mod_simu_R01WILD500_P05_M0016
mod_simu_R12 = mod_simu_R02WILD500_P05_M0016
mod_simu_R13 = mod_simu_R03WILD500_P05_M0016
mod_simu_R14 = mod_simu_R04WILD500_P05_M0016
mod_simu_R15 = mod_simu_R05WILD500_P05_M0016
mod_simu_R16 = mod_simu_R06WILD500_P05_M0016
mod_simu_R17 = mod_simu_R07WILD500_P05_M0016
mod_simu_R18 = mod_simu_R08WILD500_P05_M0016
mod_simu_R19 = mod_simu_R09WILD500_P05_M0016
mod_simu_R20 = mod_simu_R10WILD500_P05_M0016

mod_simu_R21 = mod_simu_R01WILD500_P05_M0080
mod_simu_R22 = mod_simu_R02WILD500_P05_M0080
mod_simu_R23 = mod_simu_R03WILD500_P05_M0080
mod_simu_R24 = mod_simu_R04WILD500_P05_M0080
mod_simu_R25 = mod_simu_R05WILD500_P05_M0080
mod_simu_R26 = mod_simu_R06WILD500_P05_M0080
mod_simu_R27 = mod_simu_R07WILD500_P05_M0080
mod_simu_R28 = mod_simu_R08WILD500_P05_M0080
mod_simu_R29 = mod_simu_R09WILD500_P05_M0080
mod_simu_R30 = mod_simu_R10WILD500_P05_M0080

mod_simu_R31 = mod_simu_R01WILD500_P05
mod_simu_R32 = mod_simu_R02WILD500_P05
mod_simu_R33 = mod_simu_R03WILD500_P05
mod_simu_R34 = mod_simu_R04WILD500_P05
mod_simu_R35 = mod_simu_R05WILD500_P05
mod_simu_R36 = mod_simu_R06WILD500_P05
mod_simu_R37 = mod_simu_R07WILD500_P05
mod_simu_R38 = mod_simu_R08WILD500_P05
mod_simu_R39 = mod_simu_R09WILD500_P05
mod_simu_R40 = mod_simu_R10WILD500_P05

mod_simu_R41 = mod_simu_R01WILD500_P10_M0016
mod_simu_R42 = mod_simu_R02WILD500_P10_M0016
mod_simu_R43 = mod_simu_R03WILD500_P10_M0016
mod_simu_R44 = mod_simu_R04WILD500_P10_M0016
mod_simu_R45 = mod_simu_R05WILD500_P10_M0016
mod_simu_R46 = mod_simu_R06WILD500_P10_M0016
mod_simu_R47 = mod_simu_R07WILD500_P10_M0016
mod_simu_R48 = mod_simu_R08WILD500_P10_M0016
mod_simu_R49 = mod_simu_R09WILD500_P10_M0016
mod_simu_R50 = mod_simu_R10WILD500_P10_M0016

mod_simu_R51 = mod_simu_R01WILD500_P10_M0032
mod_simu_R52 = mod_simu_R02WILD500_P10_M0032
mod_simu_R53 = mod_simu_R03WILD500_P10_M0032
mod_simu_R54 = mod_simu_R04WILD500_P10_M0032
mod_simu_R55 = mod_simu_R05WILD500_P10_M0032
mod_simu_R56 = mod_simu_R06WILD500_P10_M0032
mod_simu_R57 = mod_simu_R07WILD500_P10_M0032
mod_simu_R58 = mod_simu_R08WILD500_P10_M0032
mod_simu_R59 = mod_simu_R09WILD500_P10_M0032
mod_simu_R60 = mod_simu_R10WILD500_P10_M0032

mod_simu_R61 = mod_simu_R01WILD500_P10_M0160
mod_simu_R62 = mod_simu_R02WILD500_P10_M0160
mod_simu_R63 = mod_simu_R03WILD500_P10_M0160
mod_simu_R64 = mod_simu_R04WILD500_P10_M0160
mod_simu_R65 = mod_simu_R05WILD500_P10_M0160
mod_simu_R66 = mod_simu_R06WILD500_P10_M0160
mod_simu_R67 = mod_simu_R07WILD500_P10_M0160
mod_simu_R68 = mod_simu_R08WILD500_P10_M0160
mod_simu_R69 = mod_simu_R09WILD500_P10_M0160
mod_simu_R70 = mod_simu_R10WILD500_P10_M0160

mod_simu_R71 = mod_simu_R01WILD500_P10_M0320
mod_simu_R72 = mod_simu_R02WILD500_P10_M0320
mod_simu_R73 = mod_simu_R03WILD500_P10_M0320
mod_simu_R74 = mod_simu_R04WILD500_P10_M0320
mod_simu_R75 = mod_simu_R05WILD500_P10_M0320
mod_simu_R76 = mod_simu_R06WILD500_P10_M0320
mod_simu_R77 = mod_simu_R07WILD500_P10_M0320
mod_simu_R78 = mod_simu_R08WILD500_P10_M0320
mod_simu_R79 = mod_simu_R09WILD500_P10_M0320
mod_simu_R80 = mod_simu_R10WILD500_P10_M0320

rm(mod_simu_R01WILD500_P05_M0008,
   mod_simu_R01WILD500_P05_M0016,
   mod_simu_R01WILD500_P05_M0080,
   mod_simu_R02WILD500_P05_M0008,
   mod_simu_R02WILD500_P05_M0016,
   mod_simu_R02WILD500_P05_M0080,
   mod_simu_R03WILD500_P05_M0008,
   mod_simu_R03WILD500_P05_M0016,
   mod_simu_R03WILD500_P05_M0080,
   mod_simu_R04WILD500_P05_M0008,
   mod_simu_R04WILD500_P05_M0016,
   mod_simu_R04WILD500_P05_M0080,
   mod_simu_R05WILD500_P05_M0008,
   mod_simu_R05WILD500_P05_M0016,
   mod_simu_R05WILD500_P05_M0080,
   mod_simu_R06WILD500_P05_M0008,
   mod_simu_R06WILD500_P05_M0016,
   mod_simu_R06WILD500_P05_M0080,
   mod_simu_R07WILD500_P05_M0008,
   mod_simu_R07WILD500_P05_M0016,
   mod_simu_R07WILD500_P05_M0080,
   mod_simu_R08WILD500_P05_M0008,
   mod_simu_R08WILD500_P05_M0016,
   mod_simu_R08WILD500_P05_M0080,
   mod_simu_R09WILD500_P05_M0008,
   mod_simu_R09WILD500_P05_M0016,
   mod_simu_R09WILD500_P05_M0080,
   mod_simu_R10WILD500_P05_M0008,
   mod_simu_R10WILD500_P05_M0016,
   mod_simu_R10WILD500_P05_M0080,
   mod_simu_R01WILD500_P05,
   mod_simu_R02WILD500_P05,
   mod_simu_R03WILD500_P05,
   mod_simu_R04WILD500_P05,
   mod_simu_R05WILD500_P05,
   mod_simu_R06WILD500_P05,
   mod_simu_R07WILD500_P05,
   mod_simu_R08WILD500_P05,
   mod_simu_R09WILD500_P05,
   mod_simu_R10WILD500_P05,
   mod_simu_R01WILD500_P10_M0016,
   mod_simu_R01WILD500_P10_M0032,
   mod_simu_R01WILD500_P10_M0160,
   mod_simu_R01WILD500_P10_M0320,
   mod_simu_R02WILD500_P10_M0016,
   mod_simu_R02WILD500_P10_M0032,
   mod_simu_R02WILD500_P10_M0160,
   mod_simu_R02WILD500_P10_M0320,
   mod_simu_R03WILD500_P10_M0016,
   mod_simu_R03WILD500_P10_M0032,
   mod_simu_R03WILD500_P10_M0160,
   mod_simu_R03WILD500_P10_M0320,
   mod_simu_R04WILD500_P10_M0016,
   mod_simu_R04WILD500_P10_M0032,
   mod_simu_R04WILD500_P10_M0160,
   mod_simu_R04WILD500_P10_M0320,
   mod_simu_R05WILD500_P10_M0016,
   mod_simu_R05WILD500_P10_M0032,
   mod_simu_R05WILD500_P10_M0160,
   mod_simu_R05WILD500_P10_M0320,
   mod_simu_R06WILD500_P10_M0016,
   mod_simu_R06WILD500_P10_M0032,
   mod_simu_R06WILD500_P10_M0160,
   mod_simu_R06WILD500_P10_M0320,
   mod_simu_R07WILD500_P10_M0016,
   mod_simu_R07WILD500_P10_M0032,
   mod_simu_R07WILD500_P10_M0160,
   mod_simu_R07WILD500_P10_M0320,
   mod_simu_R08WILD500_P10_M0016,
   mod_simu_R08WILD500_P10_M0032,
   mod_simu_R08WILD500_P10_M0160,
   mod_simu_R08WILD500_P10_M0320,
   mod_simu_R09WILD500_P10_M0016,
   mod_simu_R09WILD500_P10_M0032,
   mod_simu_R09WILD500_P10_M0160,
   mod_simu_R09WILD500_P10_M0320,
   mod_simu_R10WILD500_P10_M0016,
   mod_simu_R10WILD500_P10_M0032,
   mod_simu_R10WILD500_P10_M0160,
   mod_simu_R10WILD500_P10_M0320)

Valist = list()
Vrlist = list()

for (i in 1:80){
  Valist[[i]] = get(paste("mod_simu_R",i,sep = ""))$VCV[,"animal"]
  Vrlist[[i]] = get(paste("mod_simu_R",i,sep = ""))$VCV[,"units"]
}

Vplist = Map("+", Valist,Vrlist)
h2list = Map("/",Valist,Vplist)

Vas =lapply(Valist, function(x) median(x))
Vrs =lapply(Vrlist, function(x) median(x))
Vps =lapply(Vplist, function(x) median(x))
h2s =lapply(h2list, function(x) median(x))

Vas.ci =lapply(Valist, function(x) HPDinterval(x))
Vrs.ci =lapply(Vrlist, function(x) HPDinterval(x))
Vps.ci =lapply(Vplist, function(x) HPDinterval(x))
h2s.ci =lapply(h2list, function(x) HPDinterval(x))
####stop here?
QGtable = data.frame(Va = NA,
                     Va.LCI = NA,
                     Va.UCI = NA,
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
  QGtable[i,4] <- unlist(Vrs[[i]])
  QGtable[i,5] <- unlist(Vrs.ci[[i]])[,1]
  QGtable[i,6] <- unlist(Vrs.ci[[i]])[,2]
  QGtable[i,7] <- unlist(Vps[[i]])
  QGtable[i,8] <- unlist(Vps.ci[[i]])[,1]
  QGtable[i,9] <- unlist(Vps.ci[[i]])[,2]
  QGtable[i,10] <- unlist(h2s[[i]])
  QGtable[i,11] <- unlist(h2s.ci[[i]])[,1]
  QGtable[i,12] <- unlist(h2s.ci[[i]])[,2]
}

####### plot QG params
library(tidyverse)
library(PNWColors)
pal = pnw_palette(name="Sunset",n=7,type="discrete")

QGtable$rep = paste0(rep("R", length(filenames)), rep(1:length(filenames), 1))


#### gather
QG1 =  QGtable %>% gather(varcomp, mode, c("Va","Vr","Vp","h2"))
QG2 =  QGtable %>% gather(varcomp, ci.low, c("Va.LCI","Vr.LCI","Vp.LCI","h2.LCI"))
QG3 =  QGtable %>% gather(varcomp, ci.up, c("Va.UCI","Vr.UCI","Vp.UCI","h2.UCI"))

QG = data.frame(rep = QG1$rep,
                varcomp = factor(QG1$varcomp),
                mode = QG1$mode,
                ci.low = QG2$ci.low,
                ci.up = QG3$ci.up)

QG$type = factor(rep("Wild",80))
QG$n = factor(rep(c(rep("_P05_M0008",10),rep("_P05_M0016",10),rep("_P05_M0080",10),
                    rep("_P05_M0160",10),
                    rep("_P10_M0016",10),rep("_P10_M0032",10),rep("_P10_M0160",10),
                    rep("_P10_M0320",10)),4)) # WILD

#QGfinal = read.table("QGtable.txt",header=T,stringsAsFactors = T)S
QG$type_n = factor(paste(QG$type,QG$n,sep=""))

sum.table = QG %>% 
  group_by(varcomp,type_n) %>% 
  summarize(across(everything(), list(mean)))

allplot.mode=
  sum.table %>%
  ggplot(aes(x = type_n, y = mode_1, color = varcomp))+
  geom_point(size=3) +
  geom_errorbar(aes(ymin = ci.low_1, ymax = ci.up_1,width=.02), size=1) +
  theme_bw() +
  facet_wrap(~varcomp,ncol = 2,scales="free")+
  #geom_hline(data = lines2, aes(yintercept = Z))+
  guides(color="none")+
  scale_x_discrete(limits=c("Wild_P05_M0008","Wild_P05_M0016","Wild_P05_M0080",
                            "Wild_P10_M0016","Wild_P10_M0032","Wild_P10_M0160","Wild_P10_M0320"))+
  scale_color_manual(values=pal[1:5],
                     aesthetics=c("fill","color")) +
  theme_bw() +
  ylab(expression(paste("")))+
  xlab(expression(paste("")))+
  #ggtitle(label = "Simulation study | Comparison of pedigree structure")+
  theme(#aspect.ratio=0.25,
    legend.position = "none",
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank())+
  #panel.background = element_rect(fill = "#F3F1E7", colour =NA, size = 0.5, linetype = "solid"))+
  theme(axis.text.x=element_text(angle=45, hjust=1,face = "bold",size=6))

allplot.mode+ 
  theme(aspect.ratio = 1)
#allplot.mode + geom_rect(data=QG,aes(xmin = -Inf, xmax = Inf, ymin = shadowmin, ymax = shadowmax),
#          alpha = 0.01,
#          fill = "red")


### COMPARE WITH WHAT WAS SIMULATED

simsadd = dir(pattern =".txt", full.names = TRUE)

va = {}
vp = {}
h2 = {}
va.sd = {}
vp.sd = {}
h2.sd = {}
for(i in simsadd) {
  dat = read.table(i, header = TRUE)
  dat = dat[dat$generation==10000,]
  va = c(va, mean(dat$off.q1.Va))
  va.sd = c(va.sd, sd(dat$off.q1.Va))
  vp = c(vp, mean(dat$off.q1.Vp))
  vp.sd = c(vp.sd, sd(dat$off.q1.Vp))
  h2 = c(h2, mean(dat$off.q1.Va/dat$off.q1.Vp))
  h2.sd = c(h2.sd, sd(dat$off.q1.Va/dat$off.q1.Vp))
}

###

Vas.sd =lapply(Valist, function(x) sd(x))
Vrs.sd =lapply(Vrlist, function(x) sd(x))
Vps.sd =lapply(Vplist, function(x) sd(x))
h2s.sd =lapply(h2list, function(x) sd(x))

QGtable = data.frame(Va = NA,
                     Va.sd = NA,
                     Vr = NA,
                     Vr.sd = NA,
                     Vp = NA,
                     Vp.sd = NA,
                     h2 = NA,
                     h2.sd = NA)

for (i in 1:length(Vas)) {
  QGtable[i,1] <- unlist(Vas[[i]])
  QGtable[i,2] <- unlist(Vas.sd[[i]])
  QGtable[i,3] <- unlist(Vrs[[i]])
  QGtable[i,4] <- unlist(Vrs.sd[[i]])
  QGtable[i,5] <- unlist(Vps[[i]])
  QGtable[i,6] <- unlist(Vps.sd[[i]])
  QGtable[i,7] <- unlist(h2s[[i]])
  QGtable[i,8] <- unlist(h2s.sd[[i]])
}

####### plot QG params

QGtable$rep = factor(c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10",
                       "R11","R12","R13","R14","R15","R16","R17","R18","R19","R20",
                       "R21","R22","R23","R24","R25","R26","R27","R28","R29","R30",
                       "R31","R32","R33","R34","R35","R36","R37","R38","R39","R40",
                       "R41","R42","R43","R44","R45","R46","R47","R48","R49","R50",
                       "R51","R52","R53","R54","R55","R56","R57","R58","R59","R60",
                       "R61","R62","R63","R64","R65","R66","R67","R68","R69","R70",
                       "R71","R72","R73","R74","R75","R76","R77","R78","R79","R80"))

#### gather
QG1 =  QGtable %>% gather(varcomp, mode, c("Va","Vr","Vp","h2"))
QG2 =  QGtable %>% gather(varcomp, sd, c("Va.sd","Vr.sd","Vp.sd","h2.sd"))

QG = data.frame(rep = QG1$rep,
                varcomp = factor(QG1$varcomp),
                mode = QG1$mode,
                sd = QG2$sd)

QG$type = factor(rep("Wild",80))
QG$n = factor(rep(c(rep("_P05_M0008",10),rep("_P05_M0016",10),rep("_P05_M0080",10),
                    rep("_P05_M0160",10),
                    rep("_P10_M0016",10),rep("_P10_M0032",10),rep("_P10_M0160",10),
                    rep("_P10_M0320",10)),4)) # WILD

QG$type_n = factor(paste(QG$type,QG$n,sep=""))

sum.table = QG %>% 
  group_by(varcomp,type_n) %>% 
  summarize(across(everything(), list(mean)))

names(sum.table)[5]="sd"

###
colz=c(1,2,4,5)
data.emp=sum.table[,colz]

vas=cbind(va,va.sd)
vps=cbind(vp,vp.sd)
h2s=cbind(h2,h2.sd)

data.sim=as.data.frame(rbind(vas,vps,h2s))
colnames(data.sim)=c("mode_1","sd")

data.sim$varcomp=c(rep("Va",8), rep("Vp",8),rep("h2",8))
data.sim$type_n=rep(c("Wild_P05_M0008","Wild_P05_M0016","Wild_P05_M0080","Wild_P05_M0160",
                      "Wild_P10_M0016","Wild_P10_M0032","Wild_P10_M0160","Wild_P10_M0320"),3)

toto=as.data.frame(rbind(data.emp,data.sim))
toto=toto[!toto$varcomp=="Vr",]

toto$data.type=c(rep("empirical",24),rep("simulated",24))

toto=toto[!toto$varcomp=="Vp",]

allplot.mode=
  toto %>%
  ggplot(aes(x = type_n, y = mode_1, color = data.type))+
  geom_point(size=2, alpha=0.5) +
  geom_errorbar(aes(ymin = mode_1-sd, ymax = mode_1+sd, width=.05), size=1,alpha=0.5) +
  theme_bw() +
  facet_wrap(~varcomp,ncol = 2,scales="free")+
  #geom_hline(data = lines2, aes(yintercept = Z))+
  guides(color="none")+
  scale_x_discrete(limits=c("Wild_P05_M0008","Wild_P05_M0016","Wild_P05_M0080","Wild_P05_M0160",
                            "Wild_P10_M0016","Wild_P10_M0032","Wild_P10_M0160","Wild_P10_M0320"))+
  scale_color_manual(values=c("red","blue"),
                     aesthetics=c("fill","color")) +
  theme_bw() +
  ylab(expression(paste("")))+
  xlab(expression(paste("")))+
  #ggtitle(label = "Simulation study | Comparison of pedigree structure")+
  theme(#aspect.ratio=0.25,
    legend.position = "none",
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank())+
  #panel.background = element_rect(fill = "#F3F1E7", colour =NA, size = 0.5, linetype = "solid"))+
  theme(axis.text.x=element_text(angle=45, hjust=1,face = "bold",size=6))

Fig3=allplot.mode+ 
  theme(aspect.ratio = 1)