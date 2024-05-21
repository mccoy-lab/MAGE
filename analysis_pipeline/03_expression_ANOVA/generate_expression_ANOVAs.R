# LOAD PACKAGES
library(DESeq2)
library(tidyverse)
library(nlme)
library(broom)
library(lme4)
library(broom)

# LOAD DATA
VST_expression<-read.csv("expression.vst.csv",
                         row.names=1) # MAGE.v1.0.data/global_trend_results/global_expression_trends/expression.vst.csv 

VST_metadata<-read.table("sample.metadata.MAGE.v1.0.txt",
                         header = T,
                         row.names = "sample_coriellID") # MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt

high_expression_genes<-read.table("expression_filteredGenes.MAGE.v1.0.txt.gz",
                                  header=T) %>%
  pull(ensemblID) #MAGE.v1.0.data/QTL_results/eQTL_results/expression_filteredGenes.MAGE.v1.0.txt.gz

## OBJECT PREPARATION
Expression_matrix.complete<-data.frame(VST_expression)
Expression_matrix.complete<-Expression_matrix.complete[rownames(Expression_matrix.complete) %in% high_expression_genes,]

Expression_matrix.complete.metadata<-data.frame(VST_metadata) %>% 
  select(batch,continentalGroup,population,sex)

Complete_ANOVA.ssq<-tibble(Gene=character(),
                           stage1.SSR_sex=numeric(),
                           stage1.SSR_batch=numeric(),
                           stage1.SSR_res=numeric(),
                           stage2a.SSR_continentalGroup=numeric(),
                           stage2a.SSR_res=numeric(),
                           stage2b.SSR_population=numeric(),
                           stage2b.SSR_res=numeric())
Complete_ANOVA.continentalGroup_var<-tibble(Gene=character(),
                                    subset=factor(),
                                    variance=numeric())
Complete_ANOVA.population_var<-tibble(Gene=character(),
                                  subset=factor(),
                                  variance=numeric())

n_genes<-nrow(Expression_matrix.complete)
i=1
for (gene in rownames(Expression_matrix.complete)){
  print(paste0("Processing: ",i," of ",n_genes,sep=""))
  
  #Prep gene dataframe for regression
  df<-as.data.frame(t(Expression_matrix.complete[gene,]))
  colnames(df)<-c("expression")
  df$population<-Expression_matrix.complete.metadata[rownames(df),]$population
  df$continentalGroup<-Expression_matrix.complete.metadata[rownames(df),]$continentalGroup
  df$batch<-as.factor(Expression_matrix.complete.metadata[rownames(df),]$batch)
  df$sex<-as.factor(Expression_matrix.complete.metadata[rownames(df),]$sex)
  
  #Run anova for gene
  
  #Two-step regression
  #Full array
  ssq.twostep.one<-anova(lm(expression~sex+batch,data=df))$"Sum Sq"
  names(ssq.twostep.one)<-c("sex.twostep","batch.twostep","res.twostep.first")
  df$residuals<-summary(lm(expression~sex+batch,data=df))$residuals
  
  ssq.twostep.continentalGroup<-anova(lm(residuals~continentalGroup,data=df))$"Sum Sq"
  names(ssq.twostep.continentalGroup)<-c("continentalGroup.twostep","res.twostep.continentalGroup")
  ssq.twostep.population<-anova(lm(residuals~population,data=df))$"Sum Sq"
  names(ssq.twostep.population)<-c("population.twostep","res.twostep.population")
  
  ssq<-c(ssq.twostep.one,ssq.twostep.continentalGroup,ssq.twostep.population)
  
  Complete_ANOVA.ssq<-Complete_ANOVA.ssq %>% 
    add_row(Gene=gene,
            stage1.SSR_sex=ssq["sex.twostep"],
            stage1.SSR_batch=ssq["batch.twostep"],
            stage1.SSR_res=ssq["res.twostep.first"],
            stage2a.SSR_continentalGroup=ssq["continentalGroup.twostep"],
            stage2a.SSR_res=ssq["res.twostep.continentalGroup"],
            stage2b.SSR_population=ssq["population.twostep"],
            stage2b.SSR_res=ssq["res.twostep.population"])
  i=i+1
} #All pops ANOVA
i=1
for (gene in rownames(Expression_matrix.complete)){
  print(paste0("Processing: ",i," of ",n_genes,sep=""))
  
  #Prep gene dataframe for regression
  df<-as.data.frame(t(Expression_matrix.complete[gene,]))
  colnames(df)<-c("expression")
  df$population<-Expression_matrix.complete.metadata[rownames(df),]$population
  df$continentalGroup<-Expression_matrix.complete.metadata[rownames(df),]$continentalGroup
  df$batch<-as.factor(Expression_matrix.complete.metadata[rownames(df),]$batch)
  df$sex<-as.factor(Expression_matrix.complete.metadata[rownames(df),]$sex)
  
  #Run anova for gene
  #Regress out sex and batch effects
  df$residuals<-summary(lm(expression~sex+batch,data=df))$residuals
  for (pop in c("AFR","EUR","SAS","EAS","AMR")){
    Complete_ANOVA.continentalGroup_var<-Complete_ANOVA.continentalGroup_var %>%
      add_row(Gene=gene,
              subset=pop,
              variance=(df %>% 
                          filter(continentalGroup == pop) %>% 
                          pull(residuals) %>% 
                          var()))
  }
  i=i+1
} #Partitioned ANOVAs for continentalGroup variance
i=1
for (gene in rownames(Expression_matrix.complete)){
  print(paste0("Processing: ",i," of ",n_genes,sep=""))
  
  #Prep gene dataframe for regression
  df<-as.data.frame(t(Expression_matrix.complete[gene,]))
  colnames(df)<-c("expression")
  df$population<-Expression_matrix.complete.metadata[rownames(df),]$population
  df$continentalGroup<-Expression_matrix.complete.metadata[rownames(df),]$continentalGroup
  df$batch<-as.factor(Expression_matrix.complete.metadata[rownames(df),]$batch)
  df$sex<-as.factor(Expression_matrix.complete.metadata[rownames(df),]$sex)
  
  #Run anova for gene
  #Regress out sex and batch effects
  df$residuals<-summary(lm(expression~sex+batch,data=df))$residuals
  for (pop in unique(Expression_matrix.complete.metadata$population)){
    Complete_ANOVA.population_var<-Complete_ANOVA.population_var %>%
      add_row(Gene=gene,
              subset=pop,
              variance=(df %>%
                          filter(population == pop) %>%
                          pull(residuals) %>%
                          var()))
  }
  i=i+1
} #Partitioned ANOVAs for population variance

Complete_ANOVA<-list("TwoStage_SSQ"=(Complete_ANOVA.ssq %>% 
       mutate(continentalGroupPVE=stage2a.SSR_continentalGroup/(stage2a.SSR_continentalGroup+stage2a.SSR_res),
              populationPVE=stage2b.SSR_population/(stage2b.SSR_population+stage2b.SSR_res))),
     "continentalGroup_Var"=(Complete_ANOVA.continentalGroup_var %>%
                       group_by(subset)),
     "population_Var"=(Complete_ANOVA.population_var %>%
                     group_by(subset)))

Complete_ANOVA$population_Var<-Complete_ANOVA$population_Var %>%
  left_join((Expression_matrix.complete.metadata %>% tibble() %>% select(continentalGroup,population) %>% distinct()),by=c("subset"="population")) %>%
  select(Gene,subset,continentalGroup,variance) %>%
  group_by(continentalGroup,subset)
