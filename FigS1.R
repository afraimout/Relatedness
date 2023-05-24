### Fig. S1: visualize simulated wild pop relatedness
setwd("~/../Desktop/reltedness_2k22/DRYAD/results/Simulations/MCMCglmm/")
# load representative a): 5 sub-pop & low migration
load("datafiles_for_R01WILD500_P10_M0016_addQTL.RData")

# get square GRM from dgCMatrix
GRM = as.matrix(Gadd)
GRM_square = solve(GRM)

image(GRM_square)#plot the heatmap

# load representative b): 5 sub-pop & high migration
load("datafiles_for_R01WILD500_P10_M0320_addQTL.RData")

# get square GRM from dgCMatrix
GRM = as.matrix(Gadd)
GRM_square = solve(GRM)

image(GRM_square)#plot the heatmap
