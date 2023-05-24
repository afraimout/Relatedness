### Script for Fig. S3
setwd("~/../Desktop/reltedness_2k22/DRYAD/results/Chromosome_partitioning/")
library(PNWColors); library(ggrepel)

#load chromosome partitioning results
Kpart_HR = read.table("chrom_part_HELRYT_SIZE.txt",header=T,sep="\t")
Kpart_HP = read.table("chrom_part_HELPYO_SIZE.txt",header=T,sep="\t")
Kpart_HB = read.table("chrom_part_HELBYN_SIZE.txt",header=T,sep="\t")
Kpart_HEL = read.table("chrom_part_HELSINKI_SIZE.txt",header=T,sep="\t")

HRSL = droplevels(subset(Kpart_HR, trait=="SL"))
cor.test(HRSL$size,HRSL$h2)
HRPL = droplevels(subset(Kpart_HR, trait=="PL"))
cor.test(HRPL$size,HRPL$h2)
HRBD = droplevels(subset(Kpart_HR, trait=="BD"))
cor.test(HRBD$size,HRPL$h2)

HPSL = droplevels(subset(Kpart_HP, trait=="SL"))
cor.test(HPSL$size,HPSL$h2)
HPPL = droplevels(subset(Kpart_HP, trait=="PL"))
cor.test(HPPL$size,HPPL$h2)
HPBD = droplevels(subset(Kpart_HP, trait=="BD"))
cor.test(HPBD$size,HPBD$h2)

HBSL = droplevels(subset(Kpart_HB, trait=="SL"))
cor.test(HBSL$size,HBSL$h2)
HBPL = droplevels(subset(Kpart_HB, trait=="PL"))
cor.test(HBPL$size,HBPL$h2)
HBBD = droplevels(subset(Kpart_HB, trait=="BD"))
cor.test(HBBD$size,HBBD$h2)

HELSL = droplevels(subset(Kpart_HEL, trait=="SL"))
cor.test(HELSL$size,HELSL$h2)
HELPL = droplevels(subset(Kpart_HEL, trait=="PL"))
cor.test(HELPL$size,HELPL$h2)
HELBD = droplevels(subset(Kpart_HEL, trait=="BD"))
cor.test(HELBD$size,HELBD$h2)

library(ggplot2)

theme = theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

######## QTL plots

pal <- pnw_palette("Sailboat",n=21)

HR = ggplot(data=Kpart_HR,aes(x=size, y=h2,fill=factor(LG),color=factor(LG)))+
  geom_point(size=5)+
  theme+
  facet_wrap(~trait, scales = "free_y")+
  labs(y ="Heritability", x = "")+
  geom_text_repel(aes(label=LG),size = 4,color="black",fontface="bold",max.overlaps = Inf) +
  guides(fill="none")+
  guides(color="none")+
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(#aspect.ratio=0.25,
    legend.position = "none",
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.background = element_rect(fill = "seashell", colour =NA, linewidth = 0.5, linetype = "solid"))+
  theme(axis.text.x=element_text(face = "bold"))

HP = ggplot(data=Kpart_HP,aes(x=size, y=h2,fill=factor(LG),color=factor(LG)))+
  geom_point(size=5)+
  theme+
  facet_wrap(~trait, scales = "free_y")+
  labs(y ="Heritability", x = "")+
  geom_text_repel(aes(label=LG),size = 4,color="black",fontface="bold",max.overlaps = Inf) +
  guides(fill="none")+
  guides(color="none")+
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(#aspect.ratio=0.25,
    legend.position = "none",
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.background = element_rect(fill = "seashell", colour =NA, linewidth = 0.5, linetype = "solid"))+
  theme(axis.text.x=element_text(face = "bold"))

HB = ggplot(data=Kpart_HB,aes(x=size, y=h2,fill=factor(LG),color=factor(LG)))+
  geom_point(size=5)+
  theme+
  facet_wrap(~trait, scales = "free_y")+
  labs(y ="Heritability", x = "")+
  geom_text_repel(aes(label=LG),size = 4,color="black",fontface="bold",max.overlaps = Inf) +
  guides(fill="none")+
  guides(color="none")+
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(#aspect.ratio=0.25,
    legend.position = "none",
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.background = element_rect(fill = "seashell", colour =NA, linewidth = 0.5, linetype = "solid"))+
  theme(axis.text.x=element_text(face = "bold"))


HEL = ggplot(data=Kpart_HEL,aes(x=size, y=h2,fill=factor(LG),color=factor(LG)))+
  geom_point(size=5)+
  theme+
  facet_wrap(~trait, scales = "free_y")+
  labs(y ="Heritability", x = "Chromosome length (Morgan)")+
  geom_text_repel(aes(label=LG),size = 4,color="black",fontface="bold",max.overlaps = Inf) +
  guides(fill="none")+
  guides(color="none")+
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(#aspect.ratio=0.25,
    legend.position = "none",
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.background = element_rect(fill = "seashell", colour =NA, linewidth = 0.5, linetype = "solid"))+
  theme(axis.text.x=element_text(face = "bold"))

source("../../../multiplot.R")

plot = multiplot(HR+ggtitle("HEL x RYT"),HP+ggtitle("HEL x PYO"),cols = 1)
plot2 = multiplot(HB+ggtitle("HEL x BYN"),HEL+ggtitle("HELSINKI"),cols = 1)
