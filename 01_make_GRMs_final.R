### Make GRMs
setwd("~/../Desktop/reltedness_2k22/")
library(snpReady);library(gap);library(nadiv)
source("./DRYAD/00_prepare_data_final.R") #script to format data

# Make SNP GRM for Helsinki dataset
Hel <- raw.data(X1, frame="wide", base=FALSE, 
                sweep.sample= 0.95,
                call.rate=0.90, 
                maf=0.01, imput=TRUE, imput.type = "mean")#15247 Markers

#Building the GRMs
GRMs = G.matrix(Hel$M.clean, method="VanRaden", format="wide", plot = F)

# Define additive and dominance GRMs
GA_Hel = GRMs$Ga
GD_Hel = GRMs$Gd

# Format GRMs for MCMCglmm
N=dim(HELdata)[1] #number of individuals in the dataset
i <- rep(1:N,rep(N,N))
j <- rep(1:N,N)
s <-spMatrix(N,N,i,j,as.vector(GA_Hel))
GA<-solve(s)
class(GA) <- "dgCMatrix"
rownames(GA) <- GA@Dimnames[[1]] <- with(HELdata,animal)
rownames(GA) <- GA@Dimnames[[2]] <- with(HELdata,animal)

s <-spMatrix(N,N,i,j,as.vector(GD_Hel))
GD<-solve(s)
class(GD) <- "dgCMatrix"
rownames(GD) <- GD@Dimnames[[1]] <- with(HELdata,animal)
rownames(GD) <- GD@Dimnames[[2]] <- with(HELdata,animal)

#save(GA,GD,GA_Hel,GD_Hel,HELdata, file="DataAndGRMsHelsinki.RData")