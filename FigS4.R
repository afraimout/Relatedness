### Script for Fig S4
library(MCMCglmm)
load("./DRYAD/results/Simulations/mod_simu_FSHS10_subset_multichain.RData")
combinedchains = mcmc.list(mod_list[[1]]$VCV,
                           mod_list[[2]]$VCV,
                           mod_list[[3]]$VCV)
plot(combinedchains)

