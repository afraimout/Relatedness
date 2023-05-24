### Script for Figure S2: IBD distribution
library(sequoia)
#load the full dataset
setwd("~/../Desktop/reltedness_2k22/")
load("./DRYAD/data/DataAndGRMsHelsinki.RData")

#grab the pedigree from data
ped = HELdata[,4:6]
names(ped)[1]="id"
#reassign individual names to GRM dimnames
rownames(GA_Hel)=HELdata$animal
colnames(GA_Hel)=HELdata$animal
diag(GA_Hel)=NA #remove diagonal to avoid artefactual bars for self-relatedness values
#turn the square GRM into a dataframe of pairwise relationships
G_df = data.frame(as.table(GA_Hel))[lower.tri(GA_Hel, diag = FALSE), ]
names(G_df)=c("id.A","id.B","IBD")#rename the columns to match with below dataframes

#make a dataframe with all pairwise relationships from the pedigree
PW=ComparePairs(ped, Return = "All")$Dataframe

##isolate half-sibs (HS), full-sibs (FS) and unrelated parents
#HS
HS = droplevels(PW[as.character(PW$Ped1)=="HS",])

HS.df = inner_join(G_df, HS)

hist(HS.df$IBD)#simple histogram

#FS
FS = droplevels(PW[as.character(PW$Ped1)=="FS",])

FS.df = inner_join(G_df, FS)

hist(FS.df$IBD)#simple histogram

#Unrelated parents
rownames(GA_Hel)=HELdata$id#need to rename with original ids to extract only parents 
colnames(GA_Hel)=HELdata$id#that are unrelated, otherwise all unrelated parent-offspring 
                           #will be included

GRMpar = GA_Hel[grepl("-",rownames(GA_Hel)),
                  grepl("-",colnames(GA_Hel))]
hist(GRMpar)

Gpar_df = data.frame(as.table(GRMpar))[lower.tri(GRMpar, diag = FALSE), ]

## Plot
source("multiplot.R")

wild = ggplot(Gpar_df, aes(x=Freq))+
  geom_histogram(color="darkblue", fill ="blue", binwidth = 0.01) +
  geom_vline(xintercept=0, color="red", linetype="dashed", linewidth = 1)+
  labs(y="Number of indiv. pairs", x = "") +
  scale_y_continuous(expand = c(0,0))+
  ggtitle("Wild parents")+
  geom_segment(aes(x = min(Gpar_df$Freq), y = 500, xend = min(Gpar_df$Freq), yend = min(Gpar_df$Freq)),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=min(Gpar_df$Freq), y=610, label=round(min(Gpar_df$Freq),3),size=3.5)+
  geom_segment(aes(x = max(Gpar_df$Freq), y = 500, xend = max(Gpar_df$Freq), yend = max(Gpar_df$Freq)),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=max(Gpar_df$Freq), y=610, label=round(max(Gpar_df$Freq),3),size=3.5)+
  theme_bw()+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 12, face = "bold"))

fullsibs=ggplot(FS.df, aes(x=IBD))+
  geom_histogram(color="darkblue", fill ="purple", binwidth = 0.01) +
  geom_vline(xintercept=0.5, color="red", linetype="dashed", linewidth = 1)+
  labs(y="Number of indiv. pairs", x = "") +
  scale_y_continuous(expand = c(0,0))+
  ggtitle("Full sibs")+
  geom_segment(aes(x = min(FS.df$IBD), y = 100, xend = min(FS.df$IBD), yend = min(FS.df$IBD)),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=min(FS.df$IBD), y=120, label=round(min(FS.df$IBD),3),size=3.5)+
  geom_segment(aes(x = max(FS.df$IBD), y = 100, xend = max(FS.df$IBD), yend = max(FS.df$IBD)),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=max(FS.df$IBD), y=120, label=round(max(FS.df$IBD),3),size=3.5)+
  theme_bw()+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 12, face = "bold"))

halfsibs=ggplot(HS.df, aes(x=IBD))+
  geom_histogram(color="darkblue", fill ="lightblue", binwidth = 0.01) +
  geom_vline(xintercept=0.25, color="red", linetype="dashed", linewidth = 1)+
  labs(y="Number of indiv. pairs", x = "") +
  scale_y_continuous(expand = c(0,0))+
  ggtitle("Half sibs")+
  geom_segment(aes(x = min(HS.df$IBD), y = 140, xend = min(HS.df$IBD), yend = min(HS.df$IBD)),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=min(HS.df$IBD), y=170, label=round(min(HS.df$IBD),3),size=3.5)+
  geom_segment(aes(x = max(HS.df$IBD), y = 140, xend = max(HS.df$IBD), yend = max(HS.df$IBD)),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=max(HS.df$IBD), y=170, label=round(max(HS.df$IBD),3),size=3.5)+
  theme_bw()+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 12, face = "bold"))

FigS2 = multiplot(wild, fullsibs, halfsibs, cols = 2)
