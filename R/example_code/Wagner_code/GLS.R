### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
### Glucosinolate data from the Ecotypes experiment
####### Clear workspace ########
rm(list=ls())
####### Load source file ######
source('ecotypes_source.R')

####### Load glucosinolates data #######
gs<-as.data.frame(read.delim('raw_data//Ecotypes_field_glucosinolates.txt',sep='\t',header=TRUE)) %>%
  mutate(Genotype=factor(substr(Line,1,3))) # store Genotype (get it from Line names)

####### Process raw peak areas #######
## Standardize raw peak areas by sinigrin & sample weight (in mg) ## 
## Relative response factors from Clarke et al. 2010, and C. Olson-Manning (personal communication)
amtSinigrin<-0.05 # amount of sinigrin in each sample (micromol)
# units of GS concentrations: micromol per mg dry weight

gs<-mutate(gs,
           # standardize raw peak areas by internal standard (sinigrin) and relative response factors
           std1ME=(amtSinigrin*area_1ME)/(area_Sinigrin*1), # relative response factor = 1
           std1MP=(amtSinigrin*area_1MP)/(area_Sinigrin*1), # relative response factor = 1
           std6MSOH=(amtSinigrin*area_6MSOH)/(area_Sinigrin*1), # relative response factor = 1
           std2OH1ME=(amtSinigrin*area_2OH1ME)/(area_Sinigrin*1.32),# relative response factor = 1.32
           stdI3M=(amtSinigrin*area_I3M)/(area_Sinigrin*0.28),# relative response factor = .28
           # where weight data is available, sum all standardized glucosinolates & divide by sample weight to calculate total [GS] 
           totalGS=(std1ME+std1MP+std6MSOH+std2OH1ME+stdI3M)/Weight,
           # calculate BC-ratio (don't need sample weights, which would cancel out):
           BCrat=(std1ME+std1MP+std2OH1ME)/(std1ME+std1MP+std2OH1ME+std6MSOH))

###### Exclude endogenous plants for modeling ######
gs<-filter(gs,Genotype!='END')

####### Transformation functions: arcsine square root  #######
asr<-function(value){return(asin(sqrt(value)))}

####### Visualize distributions #######

hist(gs$BCrat,40) # bimodal, as expected
hist(asr(gs$BCrat),40) # somewhat better

hist(gs$totalGS,40) #
hist(sqrt(gs$totalGS),40) # better

####### Get age of each plant #######
# Load sample metadata from file
smd<-as.data.frame(read.delim('raw_data/plant_key.txt'))
smd$Cohort<-factor(smd$Cohort)

####### Subset leaves and roots first#######
gs$Site<-factor(gs$Site)
table(gs$Type,gs$Site)
'       Jam Mah Mil Par Sil
root     70 109  16  22  89
rosette  18  60   0  20   0
'
# Separate by organ 
gsleaf<-filter(gs,Type=='rosette') 
gsroot<-filter(gs,Type=='root') 

gsleaf<-select(smd,Plant_ID,Cohort) %>% 
  merge(.,gsleaf,by='Plant_ID',all=FALSE)
gsroot<-select(smd,Plant_ID,Cohort) %>% 
  merge(.,gsroot,by='Plant_ID',all=FALSE)

# All samples were harvestsed in 2011. Note age groups (either 3 yrs or 2 yrs)
gsleaf$Age<-as.factor(ifelse(gsleaf$Cohort=='2008','3','2'))
gsroot$Age<-as.factor(ifelse(gsroot$Cohort=='2008','3','2'))

####### model BCratio, roots #######

rootBCrat<-lmer(asr(BCrat)~Site*Genotype+Site*Age+(1|Block)+(1|Line)+(1|Batch),data=gsroot)
plot(resid(rootBCrat)~fitted(rootBCrat))
qqnorm(resid(rootBCrat));qqline(resid(rootBCrat)) 
hist(resid(rootBCrat),50)
anova(rootBCrat)
'Analysis of Variance Table of type III  with  Satterthwaite 
approximation for degrees of freedom
               Sum Sq Mean Sq NumDF   DenDF F.value  Pr(>F)  
Site          2.05733 0.51433     4  77.272 2.87930 0.02804 *
Genotype      0.80559 0.20140     4  73.725 1.12745 0.35021  
Age           0.11513 0.11513     1 259.328 0.64454 0.42281  
Site:Genotype 1.33176 0.08323    16 261.723 0.46596 0.96124  
Site:Age      0.13339 0.06669     2 256.386 0.37336 0.68879  
'
r2.LMM(rootBCrat) # 0.443182

LSM.BCrat.root<-as.data.frame(lmerTest::lsmeans(rootBCrat)[1]) # save LS means
colnames(LSM.BCrat.root)<-gsub('lsmeans.table.','',colnames(LSM.BCrat.root))

rand(rootBCrat)
"Analysis of Random effects Table:
        Chi.sq Chi.DF p.value    
Block 1.14e-13      1       1    
Line  1.78e+01      1   3e-05 ***
Batch 4.48e+01      1   2e-11 ***
" 

# Write results to file:
sink(file="tables/Table_S6_glucosinolate_stats.txt")
print("Roots, BC-ratio")
anova(rootBCrat)
r2.LMM(rootBCrat) # 0.443182
rand(rootBCrat)
sink()

#######  model totalGS, roots #######

rootTotal<-lmer(sqrt(totalGS)~Site*Genotype+Site*Age+(1|Block)+(1|Line)+(1|Batch),data=gsroot)
plot(resid(rootTotal)~fitted(rootTotal)) 
hist(resid(rootTotal),50) 
qqnorm(resid(rootTotal));qqline(resid(rootTotal)) 
anova(rootTotal)
'Analysis of Variance Table of type III  with  Satterthwaite 
approximation for degrees of freedom
                Sum Sq   Mean Sq NumDF   DenDF F.value  Pr(>F)  
Site          0.044960 0.0112401     4  49.545 2.79395 0.03603 *
Genotype      0.006339 0.0015848     4  78.559 0.39393 0.81241  
Age           0.000447 0.0004466     1  40.347 0.11100 0.74073  
Site:Genotype 0.044479 0.0027800    16 263.363 0.69101 0.80229  
Site:Age      0.007018 0.0035091     2  39.696 0.87227 0.42586  
'

r2.LMM(rootTotal) # 0.3939682
LSM.Total.root<-as.data.frame(lmerTest::lsmeans(rootTotal)[1]) # save LS means
colnames(LSM.Total.root)<-gsub('lsmeans.table.','',colnames(LSM.Total.root))

rand(rootTotal)
"Analysis of Random effects Table:
      Chi.sq Chi.DF p.value  
Block   6.01      1    0.01 *
Line    5.07      1    0.02 *
Batch   5.99      1    0.01 *
"

# write results to file:
sink(file="tables/Table_S6_glucosinolate_stats.txt", append=TRUE)
print("Roots, Total [GS]")
anova(rootTotal)
r2.LMM(rootTotal) # 0.3939682
rand(rootTotal)
sink()

####### model BCratio, leaves #######
leafBCrat<-lmer(asr(BCrat)~Site*Genotype+Site*Age+(1|Block)+(1|Line)+(1|Batch),data=gsleaf)
plot(resid(leafBCrat)~fitted(leafBCrat)) # looks good
qqnorm(resid(leafBCrat));qqline(resid(leafBCrat)) # looks great
hist(resid(leafBCrat),50) # perfect
anova(leafBCrat)
'Analysis of Variance Table of type III  with  Satterthwaite 
approximation for degrees of freedom
               Sum Sq  Mean Sq NumDF  DenDF F.value    Pr(>F)    
Site          0.04410 0.022048     2  7.497  0.7789 0.4927761    
Genotype      0.94002 0.235004     4 31.926  8.3023 0.0001041 ***
Age           0.07846 0.078456     1 51.977  2.7717 0.1019601    
Site:Genotype 0.99455 0.124319     8 54.381  4.3920 0.0003855 ***
Site:Age      0.00399 0.003985     1 51.605  0.1408 0.7090325    
---  
'
r2.LMM(leafBCrat) # 0.9547019

LSM.BCrat.leaf<-as.data.frame(lmerTest::lsmeans(leafBCrat)[1]) # save LS means
colnames(LSM.BCrat.leaf)<-gsub('lsmeans.table.','',colnames(LSM.BCrat.leaf))

rand(leafBCrat)
"Analysis of Random effects Table:
      Chi.sq Chi.DF p.value    
Block  0.000      1     1.0    
Line  43.019      1   5e-11 ***
Batch  0.345      1     0.6    
"

# write results to file:
sink(file="tables/Table_S6_glucosinolate_stats.txt", append=TRUE)
print("Leaves, BC-ratio")
anova(leafBCrat)
r2.LMM(leafBCrat) # 0.9547019
rand(leafBCrat)
sink()

####### Combine all LS means into single data frame #######
LSM.all<-rbind(mutate(LSM.BCrat.leaf,Response="BCratio",Organ="leaf",Term=rownames(LSM.BCrat.leaf)),
               mutate(LSM.BCrat.root,Response="BCratio",Organ="root",Term=rownames(LSM.BCrat.root)),
               mutate(LSM.Total.root,Response="Total",Organ="root",Term=rownames(LSM.Total.root)))%>%
  plyr::rename(replace=c('Standard.Error'='SE'))
LSM.all$Term<-factor(LSM.all$Term)

####### Figure S8b: Plot: BC ratio by genotype #######
pdf(file="plots/Fig_S8b_LSM_BCratio_Genotype.pdf",width=9,height=9)
filter(LSM.all,Response=='BCratio',grepl("^Genotype  ", Term)) %>% 
ggplot(.,aes(x=Organ,y=Estimate,color=Genotype))+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), size=1.5,width=0,position=position_dodge(width=0.3))+
  geom_point(size=3,position=position_dodge(width=0.3))+
  ylab("BC-ratio")+
  scale_colour_manual(values=popPalette)+
  theme_classic()+
  theme(legend.title= element_text(size=26),legend.text=element_text(size=24,face="bold"))+
  theme(axis.title.x=element_text(size=28,face="bold"),axis.text.x=element_text(size=24,face="bold",color=c("forest green","dark grey")))+
  theme(axis.title.y=element_text(size=28,face="bold"),axis.text.y=element_text(size=24,face="bold"))+
  theme(legend.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))
dev.off()

####### Figure S8a: Plot: BC ratio by site #######
pdf(file="plots/Fig_S8a_LSM_BCratio_Site.pdf",width=9,height=9)
filter(LSM.all,Response=='BCratio',grepl("^Site  ", Term)) %>% 
  ggplot(.,aes(x=Organ,y=Estimate,color=Site))+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), size=1.5,width=0,position=position_dodge(width=0.3))+
  geom_point(size=3,position=position_dodge(width=0.3))+
  ylab("BC-ratio")+
  scale_colour_manual(values=c(sitePalette[1:2],'forest green','purple',sitePalette[3]))+
  theme_classic()+
  theme(legend.title= element_text(size=26),legend.text=element_text(size=24,face="bold"))+
  theme(axis.title.x=element_text(size=28,face="bold"),axis.text.x=element_text(size=24,face="bold",color=c("forest green","dark grey")))+
  theme(axis.title.y=element_text(size=28,face="bold"),axis.text.y=element_text(size=24,face="bold"))+
  theme(legend.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))
dev.off()

####### Figure S8d: Plot: TotalGS by site #######
pdf(file="plots/Fig_S8d_LSM_TotalGS_Site.pdf",width=9,height=9)
filter(LSM.all,Response=='Total',grepl("^Site  ", Term)) %>% 
  ggplot(.,aes(x=Site,y=Estimate,color=Site))+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), size=1.5,width=0,position=position_dodge(width=0.3))+
  geom_point(size=3,position=position_dodge(width=0.3))+
  ylab("Total root [GLS]\n")+
  scale_colour_manual(values=c(sitePalette[1:2],'forest green','purple',sitePalette[3]))+
  theme_classic()+
  theme(legend.title= element_text(size=26),legend.text=element_text(size=24,face="bold"))+
  theme(axis.title.x=element_text(size=28,face="bold"),axis.text.x=element_text(size=24,face="bold"))+
  theme(axis.title.y=element_text(size=28,face="bold"),axis.text.y=element_text(size=20,face="bold"))+
  theme(legend.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))
dev.off()

####### Figure S8c: Plot: BC ratio GxS, leaves #######
pdf(file="plots/Fig_S8c_LSM_BCratio_GxSite.pdf",width=9,height=9)
filter(LSM.all,Response=='BCratio',grepl("^Site:Genotype  ", Term),Organ=='leaf') %>% 
  ggplot(.,aes(x=Site,y=Estimate,color=Genotype))+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), size=1.5,width=0,position=position_dodge(width=0.3))+
  geom_point(size=3,position=position_dodge(width=0.3))+
  ylab("Leaf BC-ratio")+
  scale_colour_manual(values=popPalette)+
  theme_classic()+
  theme(legend.title= element_text(size=26),legend.text=element_text(size=24,face="bold"))+
  theme(axis.title.x=element_text(size=28,face="bold"),axis.text.x=element_text(size=24,face="bold",color=c(sitePalette[1:2],'red')))+
  theme(axis.title.y=element_text(size=28,face="bold"),axis.text.y=element_text(size=24,face="bold"))+
  theme(legend.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))
dev.off()

####### Save image #######
save.image(paste0("GLS/image_",date(),".RData"))

