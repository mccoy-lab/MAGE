# LOAD LIBRARIES
library(tidyverse)
library(cowplot)
library(stats)

# LOAD DATA - FIGURE 2A
expression_ANOVA<-read.csv("PVE_regressions.csv")
variance_contGroup<-read.csv("contGroup_variance.csv")

# VARIANCE EXPLAINED BY CONTINENTAL GROUP / POPULATION LABELS
ggplot(expression_ANOVA,aes(x=continentalGroupPVE)) + 
  geom_histogram(bins=1000,alpha=0.33,fill="dodgerblue1") + 
  geom_rug(aes(x=mean(continentalGroupPVE)),color="dodgerblue1",size=1,sides="b")+
  geom_text(x=0.035,y=600,label="Continental Group",color="dodgerblue1",alpha=0.33,size=4,hjust=0)+
  geom_histogram(aes(x=populationPVE),fill="violetred2",bins=1000,alpha=0.33)+
  geom_rug(aes(x=mean(populationPVE)),color="violetred2",size=1,sides="b")+
  geom_text(x=0.105,y=200,label="Population",color="violetred2",alpha=0.33,size=4,hjust=0)+
  scale_x_continuous(limits = c(0,1),breaks=seq(0,1,0.1))+
  theme_cowplot() +
  labs(x="Proportion of variance explained by population label",y="Num. genes")

# VARIANCE WITHIN POPULATION - FIGURE 2B 
ggplot(variance_contGroup,aes(x=continentalGroup,y=variance,fill=continentalGroup)) +
  geom_violin(width=0.5,color="transparent")+
  geom_boxplot(outlier.alpha = 0,fill="white",width=0.15)+
  scale_fill_brewer(palette="Set2",guide="none")+
  scale_y_log10(limits=c(1e-3,10))+
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(linetype=2,color='gray80'))+
  labs(x="Continental group",y="Sample variance within category")
  