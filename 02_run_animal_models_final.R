######## QUANTITATIVE GENOMICS OF NINE-SPINED STICKLEBACKS ########
setwd("~/../Desktop/reltedness_2k22/DRYAD/")
load("./data/DataAndGRMsHelsinki.RData")
library(MCMCglmm)

### MCMCGLMM Models ###
HELdata$sex = factor(HELdata$sex)
HELdata$dom = HELdata$id

prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002),
                                          G2=list(V=1, nu=0.002)))

model_SL <- MCMCglmm(body_length~sex,
                     random=~id+dom,
                     ginverse=list(id=GA,dom=GD),
                     data=HELdata,
                     prior=prior,
                     burnin=3000,
                     nitt=103000,
                     thin=100,
                     verbose=TRUE)
#save(model_SL, file = "model_SL_helsinki_SNP.Rdata")


model_BD <- MCMCglmm(body_depth~sex,
                     random=~id+dom,
                     ginverse=list(id=GA,dom=GD),
                     data=HELdata,
                     prior=prior,
                     burnin=3000,
                     nitt=103000,
                     thin=100,
                     verbose=TRUE)
#save(model_BD, file = "model_BD_helsinki_SNP.Rdata")

model_PL <- MCMCglmm(pelvic_length~sex,
                     random=~id+dom,
                     ginverse=list(id=GA,dom=GD),
                     data=HELdata,
                     prior=prior,
                     burnin=3000,
                     nitt=103000,
                     thin=100,
                     verbose=TRUE)
#save(model_PL, file = "model_PL_helsinki_SNP.Rdata")