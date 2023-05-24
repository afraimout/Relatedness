library(MCMCglmm)
setwd("~/../DRYAD/results/Simulations/results/")
path=getwd()
filenames = list.files(path, pattern=".RData")

for (i in 1:length(filenames)){
  load(filenames[i])
}

mod_simu_R11 = mod_simu_R01n500
mod_simu_R12 = mod_simu_R02n500
mod_simu_R13 = mod_simu_R03n500
mod_simu_R14 = mod_simu_R04n500
mod_simu_R15 = mod_simu_R05n500
mod_simu_R16 = mod_simu_R06n500
mod_simu_R17 = mod_simu_R07n500
mod_simu_R18 = mod_simu_R08n500
mod_simu_R19 = mod_simu_R09n500
mod_simu_R20 = mod_simu_R10n500

mod_simu_R21 = mod_simu_R01n1000
mod_simu_R22 = mod_simu_R02n1000
mod_simu_R23 = mod_simu_R03n1000
mod_simu_R24 = mod_simu_R04n1000
mod_simu_R25 = mod_simu_R05n1000
mod_simu_R26 = mod_simu_R06n1000
mod_simu_R27 = mod_simu_R07n1000
mod_simu_R28 = mod_simu_R08n1000
mod_simu_R29 = mod_simu_R09n1000
mod_simu_R30 = mod_simu_R10n1000

mod_simu_R31 = mod_simu_3famR01n600_SIMU3
mod_simu_R32 = mod_simu_3famR02n600_SIMU3
mod_simu_R33 = mod_simu_3famR03n600_SIMU3
mod_simu_R34 = mod_simu_3famR04n600_SIMU3
mod_simu_R35 = mod_simu_3famR05n600_SIMU3
mod_simu_R36 = mod_simu_3famR06n600_SIMU3
mod_simu_R37 = mod_simu_3famR07n600_SIMU3
mod_simu_R38 = mod_simu_3famR08n600_SIMU3
mod_simu_R39 = mod_simu_3famR09n600_SIMU3
mod_simu_R40 = mod_simu_3famR10n600_SIMU3

mod_simu_R41 = mod_simu_5famNEWR01n1000
mod_simu_R42 = mod_simu_5famNEWR02n1000
mod_simu_R43 = mod_simu_5famNEWR03n1000
mod_simu_R44 = mod_simu_5famNEWR04n1000
mod_simu_R45 = mod_simu_5famNEWR05n1000
mod_simu_R46 = mod_simu_5famNEWR06n1000
mod_simu_R47 = mod_simu_5famNEWR07n1000
mod_simu_R48 = mod_simu_5famNEWR08n1000
mod_simu_R49 = mod_simu_5famNEWR09n1000
mod_simu_R50 = mod_simu_5famNEWR10n1000

mod_simu_R51 = mod_list[[1]]
mod_simu_R52 = mod_list[[2]]
mod_simu_R53 = mod_list[[3]]
mod_simu_R54 = mod_list[[4]]
mod_simu_R55 = mod_list[[5]]
mod_simu_R56 = mod_list[[6]]
mod_simu_R57 = mod_list[[7]]
mod_simu_R58 = mod_list[[8]]
mod_simu_R59 = mod_list[[9]]
mod_simu_R60 = mod_list[[10]]

mod_simu_R61 = mod_simu_R01HSFS03
mod_simu_R62 = mod_simu_R02HSFS03
mod_simu_R63 = mod_simu_R03HSFS03
mod_simu_R64 = mod_simu_R04HSFS03
mod_simu_R65 = mod_simu_R05HSFS03
mod_simu_R66 = mod_simu_R06HSFS03
mod_simu_R67 = mod_simu_R07HSFS03
mod_simu_R68 = mod_simu_R08HSFS03
mod_simu_R69 = mod_simu_R09HSFS03
mod_simu_R70 = mod_simu_R10HSFS03

mod_simu_R71 = mod_simu_HSFS05_R01
mod_simu_R72 = mod_simu_HSFS05_R02
mod_simu_R73 = mod_simu_HSFS05_R03
mod_simu_R74 = mod_simu_HSFS05_R04
mod_simu_R75 = mod_simu_HSFS05_R05
mod_simu_R76 = mod_simu_HSFS05_R06
mod_simu_R77 = mod_simu_HSFS05_R07
mod_simu_R78 = mod_simu_HSFS05_R08
mod_simu_R79 = mod_simu_HSFS05_R09
mod_simu_R80 = mod_simu_HSFS05_R10

mod_simu_R81 = mod_simu_R01HSFS10
mod_simu_R82 = mod_simu_R02HSFS10
mod_simu_R83 = mod_simu_R03HSFS10
mod_simu_R84 = mod_simu_R04HSFS10
mod_simu_R85 = mod_simu_R05HSFS10
mod_simu_R86 = mod_simu_R06HSFS10
mod_simu_R87 = mod_simu_R07HSFS10
mod_simu_R88 = mod_simu_R08HSFS10
mod_simu_R89 = mod_simu_R09HSFS10
mod_simu_R90 = mod_simu_R10HSFS10

mod_simu_R91 = mod_simu_R01WILD500
mod_simu_R92 = mod_simu_R02WILD500
mod_simu_R93 = mod_simu_R03WILD500
mod_simu_R94 = mod_simu_R04WILD500
mod_simu_R95 = mod_simu_R05WILD500
mod_simu_R96 = mod_simu_R06WILD500
mod_simu_R97 = mod_simu_R07WILD500
mod_simu_R98 = mod_simu_R08WILD500
mod_simu_R99 = mod_simu_R09WILD500
mod_simu_R100 = mod_simu_R10WILD500

mod_simu_R101 = mod_simu_R01WILD1000
mod_simu_R102 = mod_simu_R02WILD1000
mod_simu_R103 = mod_simu_R03WILD1000
mod_simu_R104 = mod_simu_R04WILD1000
mod_simu_R105 = mod_simu_R05WILD1000
mod_simu_R106 = mod_simu_R06WILD1000
mod_simu_R107 = mod_simu_R07WILD1000
mod_simu_R108 = mod_simu_R08WILD1000
mod_simu_R109 = mod_simu_R09WILD1000
mod_simu_R110 = mod_simu_R10WILD1000

rm(mod_simu_R01n500,
   mod_simu_R02n500,
   mod_simu_R03n500,
   mod_simu_R04n500,
   mod_simu_R05n500,
   mod_simu_R06n500,
   mod_simu_R07n500,
   mod_simu_R08n500,
   mod_simu_R09n500,
   mod_simu_R10n500,
   mod_simu_R01n1000,
   mod_simu_R02n1000,
   mod_simu_R03n1000,
   mod_simu_R04n1000,
   mod_simu_R05n1000,
   mod_simu_R06n1000,
   mod_simu_R07n1000,
   mod_simu_R08n1000,
   mod_simu_R09n1000,
   mod_simu_R10n1000,
   mod_simu_R01WILD500,
   mod_simu_R02WILD500,
   mod_simu_R03WILD500,
   mod_simu_R04WILD500,
   mod_simu_R05WILD500,
   mod_simu_R06WILD500,
   mod_simu_R07WILD500,
   mod_simu_R08WILD500,
   mod_simu_R09WILD500,
   mod_simu_R10WILD500,
   mod_simu_R01WILD1000,
   mod_simu_R02WILD1000,
   mod_simu_R03WILD1000,
   mod_simu_R04WILD1000,
   mod_simu_R05WILD1000,
   mod_simu_R06WILD1000,
   mod_simu_R07WILD1000,
   mod_simu_R08WILD1000,
   mod_simu_R09WILD1000,
   mod_simu_R10WILD1000,
   mod_simu_3famR01n600_SIMU3,
   mod_simu_3famR02n600_SIMU3,
   mod_simu_3famR03n600_SIMU3,
   mod_simu_3famR04n600_SIMU3,
   mod_simu_3famR05n600_SIMU3,
   mod_simu_3famR06n600_SIMU3,
   mod_simu_3famR07n600_SIMU3,
   mod_simu_3famR08n600_SIMU3,
   mod_simu_3famR09n600_SIMU3,
   mod_simu_3famR10n600_SIMU3,
   mod_simu_5famNEWR01n1000,
   mod_simu_5famNEWR02n1000,
   mod_simu_5famNEWR03n1000,
   mod_simu_5famNEWR04n1000,
   mod_simu_5famNEWR05n1000,
   mod_simu_5famNEWR06n1000,
   mod_simu_5famNEWR07n1000,
   mod_simu_5famNEWR08n1000,
   mod_simu_5famNEWR09n1000,
   mod_simu_5famNEWR10n1000,
   mod_simu_HSFS05_R01,
   mod_simu_HSFS05_R02,
   mod_simu_HSFS05_R03,
   mod_simu_HSFS05_R04,
   mod_simu_HSFS05_R05,
   mod_simu_HSFS05_R06,
   mod_simu_HSFS05_R07,
   mod_simu_HSFS05_R08,
   mod_simu_HSFS05_R09,
   mod_simu_HSFS05_R10,
   mod_simu_R01HSFS10,
   mod_simu_R02HSFS10,
   mod_simu_R03HSFS10,
   mod_simu_R04HSFS10,
   mod_simu_R05HSFS10,
   mod_simu_R06HSFS10,
   mod_simu_R07HSFS10,
   mod_simu_R08HSFS10,
   mod_simu_R09HSFS10,
   mod_simu_R10HSFS10,
   mod_list)

Valist = list()
Vrlist = list()

for (i in 1:110){
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

### Compare to simulated datasets

simsadd = dir(path = "../quanti/", pat=".txt", full.names = TRUE)

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
                       "R71","R72","R73","R74","R75","R76","R77","R78","R79","R80",
                       "R81","R82","R83","R84","R85","R86","R87","R88","R89","R90",
                       "R91","R92","R93","R94","R95","R96","R97","R98","R99","R100",
                       "R101","R102","R103","R104","R105","R106","R107","R108","R109","R110"))

#### gather
QG1 =  QGtable %>% gather(varcomp, mode, c("Va","Vr","Vp","h2"))
QG2 =  QGtable %>% gather(varcomp, sd, c("Va.sd","Vr.sd","Vp.sd","h2.sd"))

QG = data.frame(rep = QG1$rep,
                varcomp = factor(QG1$varcomp),
                mode = QG1$mode,
                sd = QG2$sd)

QG$type = factor(rep(c(rep("F1FS",60),rep("F1HS",30),rep("Wild",20))))
QG$n = factor(rep(c(rep(200,10),rep(500,10),rep(1000,10),rep("3x200",10),rep("5x200",10),rep("50x20",10), # FULLSIB
                    rep("150x6",10),rep("50x10",10),rep("50x20",10), # HALFSIB
                    rep(500,10),rep(1000,10)),4)) # WILD

QG$type_n = factor(paste(QG$type,QG$n,sep="_"))

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

data.sim$varcomp=c(rep("Va",11), rep("Vp",11),rep("h2",11))
data.sim$type_n=rep(c("F1FS_200","F1FS_500","F1FS_1000","F1FS_3x200","F1FS_5x200","F1FS_50x20",
                      "F1HS_150x6","F1HS_50x10","F1HS_50x20",
                      "Wild_500","Wild_1000"),3)

toto=as.data.frame(rbind(data.emp,data.sim))
toto=toto[!toto$varcomp=="Vr",]

toto$data.type=c(rep("empirical",33),rep("simulated",33))

toto=toto[!toto$varcomp=="Vp",]

allplot.mode=
  toto %>%
  ggplot(aes(x = type_n, y = mode_1, color = data.type))+
  geom_point(size=2, alpha=0.5) +
  geom_errorbar(aes(ymin = mode_1-sd, ymax = mode_1+sd, width=.05), size=1,alpha=0.5) +
  theme_bw() +
  facet_wrap(~varcomp,ncol = 2,scales="free")+
  guides(color="none")+
  scale_x_discrete(limits=c("F1FS_200","F1FS_500","F1FS_1000","F1FS_3x200","F1FS_5x200","F1FS_50x20",
                            "F1HS_50x10","F1HS_50x20","F1HS_150x6",
                            "Wild_500","Wild_1000"))+
  scale_color_manual(values=c("red","blue"),
                     aesthetics=c("fill","color")) +
  theme_bw() +
  ylab(expression(paste("")))+
  xlab(expression(paste("")))+
  theme(legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank())+
  theme(axis.text.x=element_text(angle=45, hjust=1,face = "bold",size=6))

allplot.mode+ 
  theme(aspect.ratio = 1)



