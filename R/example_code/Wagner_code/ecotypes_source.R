# Source file for Boechera_Ecotypes_Microbiome analyses

####### Color palettes #######
sitePalette <- c("#000000", "#E69F00", "#56B4E9")
popPalette <- c("#CC79A7", 
                "goldenrod", "#0072B2","#009E73","grey","#F0E442","lavender",
                "aquamarine3","darkorchid","cadetblue3","#D55E00")
agePalette<-c("lightskyblue","blue","black","red")####### Load packages #######

####### Load packages #######
library(tidyr)
library(maps)
library(mapdata)
library(ape)
library(Biostrings)
library(phyloseq)
library(genefilter)
library(plyr)
library(dplyr)
library(doParallel)
library(vegan)
library(ggplot2)
library(grid)
library(scales)
library(lme4)
library(DESeq2)
library(lmerTest)
library(reshape2)
library(BiocParallel)
library(VennDiagram)

####### Register parallel backend #######
registerDoParallel(cores=4)

####### Define zero-tolerant function for calculating geometric means (credit: P.J. McMurdie) #######
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
####### Function: r2.LMM() => estimate R^2 for a model as correlation between fitted and observed values #######
r2.LMM <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}
####### Function: LMMstats() => *with log Obs and MiSeq run as covariates*  #####
LMMstats<-function(data,resp,reml=TRUE){
  df<-data
  stats<-data.frame('Response'=character(0),'Term'=character(0), 'F or ChiSq'= numeric(0), 'df' = integer(0), 'ddf'=numeric(0),'P'=numeric(0))
  model<-lmer(df[[resp]]~Genotype*Site + Age*Site + Harvested*Site +(1|Site:Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=df,REML=reml)
  fixed<-anova(model,type=3,ddf='Satterthwaite')
  print(fixed)
  for (i in 1:length(rownames(fixed))) { # loop through fixed effects
    newrow<-data.frame('Response'=resp,'Term'=rownames(fixed)[i],'F or ChiSq'= fixed[i,5], 'df' = fixed[i,3], 'ddf'=fixed[i,4],'P'=fixed[i,6])
    stats<-rbind(stats,newrow)
  }
  ### now do random effects ### 
  BlkLRT<-rand(model)
  Blk<-data.frame('Response'=resp,'Term'='Block','F.or.ChiSq'=BlkLRT$rand.table$Chi.sq[1],"df"=BlkLRT$rand.table$Chi.DF[1],'ddf'=Inf,'P'=BlkLRT$rand.table$p.value[1])
  Line<-data.frame('Response'=resp,'Term'='Line','F.or.ChiSq'=BlkLRT$rand.table$Chi.sq[2],"df"=BlkLRT$rand.table$Chi.DF[2],'ddf'=Inf,'P'=BlkLRT$rand.table$p.value[2])
  Plate<-data.frame('Response'=resp,'Term'='Plate','F.or.ChiSq'=BlkLRT$rand.table$Chi.sq[3],"df"=BlkLRT$rand.table$Chi.DF[3],'ddf'=Inf,'P'=BlkLRT$rand.table$p.value[3])
  stats<-rbind(stats,Blk,Line,Plate)
  return(stats)
}

