#!/usr/bin/Rscript

### Script to run random-effects models for variance-stabilizing-transformed abundance of a single taxon
### Working directory should include 1) this script, and 2) "indata.RData", which is a data frame (named 'thisdf') that contains M columns of sample metadata (including linear predictors), and a Yth column with count data for one taxon.
### One taxon per working directory.
### A job array batch script will submit this program to the queue FROM EACH working directory (working dirs have integer names)

# Clear workspace 
rm(list=ls())

library(lme4)
library(lsmeans)

sessionInfo()

# This should be the same as in populate.R :
# How many columns of metadata are there? default = 19
M<-28 # first M columns of the data frame are metadata
Y<-M+1  # Taxon counts will be stored in column number Y

load('indata.RData',.GlobalEnv) # load data frame

print("Loading data...")
thisdf<-thisdata # store copy of thisdf
thistaxon<-colnames(thisdf)[Y] # store name of current taxon
save(thistaxon,file="taxon_name.RData") # for easy import

colnames(thisdf)
print(thistaxon)

print("Defining functions...")
####### lmer: LMM for variance stabilizing transformed count data #######
options(contrasts = c("contr.sum","contr.poly")) # contrasts sum to 0
lmer.model<-try(lmer((thisdf[[thistaxon]])~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=thisdf,REML=TRUE))

if (class(lmer.model)=="try-error"){
  sink("../model_errors.txt",append=TRUE)
  print(thistaxon)
  print("\n")
  sink()
}

# print diagnostic plots
pdf(file=paste(thistaxon,"--LMM_residual_plot.pdf",sep=""))
plot(residuals(lmer.model)~fitted(lmer.model),col="dark grey") # plot residuals against fitted values
abline(h=0) # add horizontal line y=0
dev.off()
pdf(file=paste(thistaxon,"--LMM_QQnorm_plot.pdf",sep=""))
qqnorm(residuals(lmer.model));qqline(residuals(lmer.model),distribution=qnorm)
dev.off()

print ("Diagnostic plots have been written to file")

####### Heritability estimates using variance components #######
# Fit random-effects model for each site
lmm.rand.Jam<-try(lmer((subset(thisdf,Site=='Jam')[[thistaxon]])~(1|Genotype)+(1|Age)+(1|Harvested)+(1|Block)+(1|Genotype:Line)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Jam'),REML=TRUE))
lmm.rand.Mah<-try(lmer((subset(thisdf,Site=='Mah')[[thistaxon]])~(1|Genotype)+(1|Age)+(1|Harvested)+(1|Block)+(1|Genotype:Line)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Mah'),REML=TRUE))
lmm.rand.Sil<-try(lmer((subset(thisdf,Site=='Sil')[[thistaxon]])~(1|Genotype)+(1|Age)+(1|Harvested)+(1|Block)+(1|Genotype:Line)+(1|newPlate)+logObs,
                  data=subset(thisdf,Site=='Sil'),REML=TRUE))
lmer.rand.all<-try(lmer(thisdf[[thistaxon]]~(1|Site)+(1|Genotype)+(1|Age)+(1|Harvested)+(1|Site:Genotype)+(1|Site:Age)+(1|Site:Harvested)+(1|Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=thisdf,REML=TRUE))

### Calculate percent variance due to genetic effects in each site 
# add % variance for Genotype and Line, divide by total 
h2.vcov.Jam<-subset(as.data.frame(VarCorr(lmm.rand.Jam)),grp=='Genotype')$vcov+
  subset(as.data.frame(VarCorr(lmm.rand.Jam)),grp=='Genotype:Line')$vcov
h2.vcov.Jam<-h2.vcov.Jam/sum(as.data.frame(VarCorr(lmm.rand.Jam))$vcov) # what proportion of the total variance?

h2.vcov.Mah<-subset(as.data.frame(VarCorr(lmm.rand.Mah)),grp=='Genotype')$vcov+
  subset(as.data.frame(VarCorr(lmm.rand.Mah)),grp=='Genotype:Line')$vcov
h2.vcov.Mah<-h2.vcov.Mah/sum(as.data.frame(VarCorr(lmm.rand.Mah))$vcov) # what proportion of the total variance?

h2.vcov.Sil<-subset(as.data.frame(VarCorr(lmm.rand.Sil)),grp=='Genotype')$vcov+
  subset(as.data.frame(VarCorr(lmm.rand.Sil)),grp=='Genotype:Line')$vcov
h2.vcov.Sil<-h2.vcov.Sil/sum(as.data.frame(VarCorr(lmm.rand.Sil))$vcov) # what proportion of the total variance?

h2.vcov.all<-subset(as.data.frame(VarCorr(lmer.rand.all)),grp=='Genotype')$vcov+
  subset(as.data.frame(VarCorr(lmer.rand.all)),grp=='Genotype:Line')$vcov+
  subset(as.data.frame(VarCorr(lmer.rand.all)),grp=='Site:Genotype')$vcov
h2.vcov.all<-h2.vcov.all/sum(as.data.frame(VarCorr(lmer.rand.all))$vcov) # what proportion of the total variance?

# initialize data frame
h2.vcov<-data.frame("Taxon"=character(),"Site"=character(),"h2_vcov"=numeric())

h2.vcov<-rbind(h2.vcov,
          data.frame("Taxon"=thistaxon,"Site"="Jam","h2_vcov"=h2.vcov.Jam),
          data.frame("Taxon"=thistaxon,"Site"="Mah","h2_vcov"=h2.vcov.Mah),
          data.frame("Taxon"=thistaxon,"Site"="Sil","h2_vcov"=h2.vcov.Sil),
          data.frame("Taxon"=thistaxon,"Site"="all","h2_vcov"=h2.vcov.all))

write.table(h2.vcov,file=paste(thistaxon,"--heritability_vcov.txt",sep=""),sep='\t',row.names=FALSE,col.names=TRUE) # save heritability estimates

print(" Heritability estimates have been written to file")
