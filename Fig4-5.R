### Script for Figure 4
setwd("~/../Desktop/reltedness_2k22/DRYAD/results/Chromosome_partitioning/")
library(ggplot2)

#load the chromosome partitioning results
Kpart_HR = read.table("chrom_part_HELRYT.txt",header=T,sep="\t")
Kpart_HR = read.table("chrom_part_HELRYT.txt",header=T,sep="\t")
Kpart_HR = read.table("chrom_part_HELRYT.txt",header=T,sep="\t")
Kpart_all = read.table("chrom_part_all.txt",header=T,sep="\t")
Kpart_HEL = read.table("chrom_part_HELSINKI.txt",header=T,sep="\t")

theme = theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

######## QTL plots

F2plot = ggplot(data=Kpart_all,aes(x=factor(LG), y=h2,fill=cross))+
  geom_bar(stat = "identity")+
  theme+
  facet_wrap(~cross+trait, scales = "free_y")+
  labs(y ="Heritability", x = "Chromosome")+
  guides(fill="none")

######## HELSINKI plots

Helsinkiplot = ggplot(data=Kpart_HEL,aes(x=factor(LG), y=h2,fill=trait))+
  geom_bar(stat = "identity")+
  theme+
  facet_wrap(trait~., scales = "free_y", ncol = 1)+
  labs(y ="Heritability", x = "Chromosome")+
  guides(fill="none")
