### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
##### Linear mixed models of alpha and beta diversity
####### Clear workspace ########
rm(list=ls())

####### Load source file #######
source('ecotypes_source.R')

####### sessionInfo #######
sessionInfo()
"R version 3.2.3 (2015-12-10)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats4    stats     graphics  grDevices utils    
[8] datasets  methods   base     

other attached packages:
[1] VennDiagram_1.6.17         futile.logger_1.4.1       
[3] BiocParallel_1.4.3         reshape2_1.4.1            
[5] lmerTest_2.0-30            lme4_1.1-11               
[7] Matrix_1.2-4               scales_0.4.0              
[9] ggplot2_2.1.0              vegan_2.3-4               
[11] lattice_0.20-33            permute_0.9-0             
[13] doParallel_1.0.10          iterators_1.0.8           
[15] foreach_1.4.3              dplyr_0.4.3               
[17] plyr_1.8.3                 genefilter_1.52.1         
[19] Biostrings_2.38.4          XVector_0.10.0            
[21] ape_3.4                    mapdata_2.2-6             
[23] maps_3.1.0                 DESeq2_1.10.1             
[25] RcppArmadillo_0.6.600.4.0  Rcpp_0.12.3               
[27] SummarizedExperiment_1.0.2 Biobase_2.30.0            
[29] GenomicRanges_1.22.4       GenomeInfoDb_1.6.3        
[31] IRanges_2.4.8              S4Vectors_0.8.11          
[33] BiocGenerics_0.16.1        phyloseq_1.14.0           

loaded via a namespace (and not attached):
[1] splines_3.2.3        Formula_1.2-1        assertthat_0.1      
[4] latticeExtra_0.6-28  RSQLite_1.0.0        digest_0.6.9        
[7] chron_2.3-47         RColorBrewer_1.1-2   minqa_1.2.4         
[10] colorspace_1.2-6     XML_3.98-1.4         zlibbioc_1.16.0     
[13] xtable_1.8-2         annotate_1.48.0      mgcv_1.8-12         
[16] lazyeval_0.1.10      nnet_7.3-12          survival_2.38-3     
[19] RJSONIO_1.3-0        magrittr_1.5         nlme_3.1-126        
[22] MASS_7.3-45          foreign_0.8-66       tools_3.2.3         
[25] data.table_1.9.6     stringr_1.0.0        munsell_0.4.3       
[28] locfit_1.5-9.1       cluster_2.0.3        AnnotationDbi_1.32.3
[31] lambda.r_1.1.7       compiler_3.2.3       ade4_1.7-4          
[34] nloptr_1.0.4         biom_0.3.12          igraph_1.0.1        
[37] labeling_0.3         gtable_0.2.0         codetools_0.2-14    
[40] multtest_2.26.0      DBI_0.3.1            R6_2.1.2            
[43] gridExtra_2.2.1      Hmisc_3.17-2         futile.options_1.0.0
[46] stringi_1.0-1        geneplotter_1.48.0   rpart_4.1-10        
[49] acepack_1.3-3.3   "

####### Load sample metadata with alpha and beta diversity scores #######
load('intermediate_data/smd_leaf3_withEndog.RData')
load('intermediate_data/smd_root3_withEndog.RData')
####### Remove endogenous plants from dataset #######
leaf3.smd<-filter(leaf3.smd.withEndog,Age!='endog')
root3.smd<-filter(root3.smd.withEndog,Age!='endog')

####### Check residuals and R^2 for all alpha diversity models and save LS means for plotting ########
leaf3.sha<-lmer(Shannon~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=leaf3.smd,REML=TRUE)
plot(residuals(leaf3.sha)~fitted(leaf3.sha)) 
qqnorm(residuals(leaf3.sha)); qqline(residuals(leaf3.sha))
LSM.sha.leaf3<-as.data.frame(lmerTest::lsmeans(leaf3.sha)[1])
colnames(LSM.sha.leaf3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(leaf3.sha) # 0.5629
rm(leaf3.sha)

leaf3.chao1<-lmer((Chao1)~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=leaf3.smd,REML=TRUE)
plot(residuals(leaf3.chao1)~fitted(leaf3.chao1))
qqnorm(residuals(leaf3.chao1)); qqline(residuals(leaf3.chao1))
LSM.chao1.leaf3<-as.data.frame(lmerTest::lsmeans(leaf3.chao1)[1])
colnames(LSM.chao1.leaf3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(leaf3.chao1) # 0.6330
rm(leaf3.chao1)

root3.sha<-lmer(Shannon~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=root3.smd,REML=TRUE)
plot(residuals(root3.sha)~fitted(root3.sha)) 
qqnorm(residuals(root3.sha)); qqline(residuals(root3.sha)) 
LSM.sha.root3<-as.data.frame(lmerTest::lsmeans(root3.sha)[1])
colnames(LSM.sha.root3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(root3.sha) # 0.6825
rm(root3.sha)

root3.chao1<-lmer(Chao1~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=root3.smd,REML=TRUE)
plot(residuals(root3.chao1)~fitted(root3.chao1)) 
qqnorm(residuals(root3.chao1)); qqline(residuals(root3.chao1))
LSM.chao1.root3<-as.data.frame(lmerTest::lsmeans(root3.chao1)[1])
colnames(LSM.chao1.root3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(root3.chao1) # 0.7184
rm(root3.chao1)

####### Table 1: Alpha diversity results #######
# Use Shannon entropy estimate as response variable because residuals are terrible when using e^Shannon 
stats.adiv.leaf3<-LMMstats(leaf3.smd,resp="Shannon",reml=TRUE) %>%
  rbind(.,LMMstats(leaf3.smd,resp="Chao1",reml=TRUE)) %>%
  Pcorrect(method="holm") 

stats.adiv.root3<-LMMstats(root3.smd,resp="Shannon",reml=TRUE) %>%
  rbind(.,LMMstats(root3.smd,resp="Chao1",reml=TRUE)) %>%
  Pcorrect(method="holm") 

### Write table to file, for further formatting in word processor
Table1<-data.frame(row.names=levels(stats.adiv.leaf3$Term))
Table1$Term<-rownames(Table1)
for (h in 1:2){
  organ<-c('leaf','root')[h]
  stats<-get(c("stats.adiv.leaf3","stats.adiv.root3")[h])
  for (i in 1:length(levels(stats.adiv.leaf3$Response))){
    resp<-levels(stats$Response)[i]
    for (j in 1:length(levels(stats$Term))){
      term<-levels(stats$Term)[j]
      subtable1<-subset(stats,Response==resp&Term==term)
      teststat<-ifelse(term%in%c('Block','Line','Plate'),
                       paste0("ChiSq",subtable1$df,"=",format(subtable1$F.or.ChiSq,digits=3)),
                       paste0("F",subtable1$df,",",ceiling(subtable1$ddf),"=",format(subtable1$F.or.ChiSq,digits=3)))
      p<-ifelse(subtable1$P_corrected==0,"P<3e-16",paste0("P=",format(subtable1$P_corrected,digits=2)))
      Table1[term,paste0(organ,"__",resp)]<-paste0(teststat,";",p)
    }
  }
}
write.table(Table1,file="tables/Table1_alphadiv_stats.txt",sep='\t',row.names=FALSE,col.names=TRUE)

####### Check residuals and R^2 for weighted UniFrac models and save LS means for plotting #######
# weighted UniFrac PCo1: leaves
leaf3.wUF1<-lmer(wUF1~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=leaf3.smd,REML=TRUE)
plot(residuals(leaf3.wUF1)~fitted(leaf3.wUF1)) 
qqnorm(residuals(leaf3.wUF1)); qqline(residuals(leaf3.wUF1)) 
LSM.wUF1.leaf3<-as.data.frame(lmerTest::lsmeans(leaf3.wUF1)[1])
colnames(LSM.wUF1.leaf3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(leaf3.wUF1) # 0.842
rm(leaf3.wUF1)

# weighted UniFrac PCo2: leaves
leaf3.wUF2<-lmer(wUF2~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=leaf3.smd,REML=TRUE)
plot(residuals(leaf3.wUF2)~fitted(leaf3.wUF2))
qqnorm(residuals(leaf3.wUF2)); qqline(residuals(leaf3.wUF2)) 
LSM.wUF2.leaf3<-as.data.frame(lmerTest::lsmeans(leaf3.wUF2)[1])
colnames(LSM.wUF2.leaf3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(leaf3.wUF2) # 0.630
rm(leaf3.wUF2)

# weighted UniFrac PCo3: leaves
leaf3.wUF3<-lmer(wUF3~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=leaf3.smd,REML=TRUE)
plot(residuals(leaf3.wUF3)~fitted(leaf3.wUF3))
qqnorm(residuals(leaf3.wUF3)); qqline(residuals(leaf3.wUF3)) 
LSM.wUF3.leaf3<-as.data.frame(lmerTest::lsmeans(leaf3.wUF3)[1])
colnames(LSM.wUF3.leaf3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(leaf3.wUF3) # 0.636
rm(leaf3.wUF3)

# weighted UniFrac PCo1: roots
root3.wUF1<-lmer(wUF1~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=root3.smd,REML=TRUE)
plot(residuals(root3.wUF1)~fitted(root3.wUF1))
qqnorm(residuals(root3.wUF1)); qqline(residuals(root3.wUF1)) 
LSM.wUF1.root3<-as.data.frame(lmerTest::lsmeans(root3.wUF1)[1])
colnames(LSM.wUF1.root3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(root3.wUF1) # 0.894
rm(root3.wUF1)

# weighted UniFrac PCo2: roots
root3.wUF2<-lmer(wUF2~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=root3.smd,REML=TRUE)
plot(residuals(root3.wUF2)~fitted(root3.wUF2))
qqnorm(residuals(root3.wUF2)); qqline(residuals(root3.wUF2)) 
LSM.wUF2.root3<-as.data.frame(lmerTest::lsmeans(root3.wUF2)[1])
colnames(LSM.wUF2.root3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(root3.wUF2) # 0.765
rm(root3.wUF2)

# weighted UniFrac PCo3: roots
root3.wUF3<-lmer(wUF3~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=root3.smd,REML=TRUE)
plot(residuals(root3.wUF3)~fitted(root3.wUF3))
qqnorm(residuals(root3.wUF3)); qqline(residuals(root3.wUF3)) 
LSM.wUF3.root3<-as.data.frame(lmerTest::lsmeans(root3.wUF3)[1])
colnames(LSM.wUF3.root3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(root3.wUF3) # 0.441
rm(root3.wUF3)


####### Table 2: weighted UniFrac results #######
# use "LMMstats" function to control for MiSeq run and sequencing depth (logObs)
stats.wUF.leaf3<-LMMstats(leaf3.smd,resp="wUF1",reml=TRUE) %>%
  rbind(.,LMMstats(leaf3.smd,resp="wUF2",reml=TRUE)) %>%
  rbind(.,LMMstats(leaf3.smd,resp="wUF3",reml=TRUE)) %>%
  Pcorrect(method='holm')

stats.wUF.root3<-LMMstats(root3.smd,resp="wUF1",reml=TRUE) %>%
  rbind(.,LMMstats(root3.smd,resp="wUF2",reml=TRUE)) %>%
  rbind(.,LMMstats(root3.smd,resp="wUF3",reml=TRUE)) %>%
  Pcorrect(method='holm')

Table2<-data.frame(row.names=levels(stats.wUF.leaf3$Term))
Table2$Term<-rownames(Table2)
for (h in 1:2){
  organ<-c('leaf','root')[h]
  stats<-get(c("stats.wUF.leaf3","stats.wUF.root3")[h])
  for (i in 1:length(levels(stats$Response))){
    resp<-levels(stats$Response)[i]
    for (j in 1:length(levels(stats$Term))){
      term<-levels(stats$Term)[j]
      subTable2<-subset(stats,Response==resp&Term==term)
      teststat<-ifelse(term%in%c('Block','Line','Plate'),
                       paste0("ChiSq",subTable2$df,"=",format(subTable2$F.or.ChiSq,digits=3)),
                       paste0("F",subTable2$df,",",ceiling(subTable2$ddf),"=",format(subTable2$F.or.ChiSq,digits=3)))
      p<-ifelse(subTable2$P_corrected==0,"P<3e-16",paste0("P=",format(subTable2$P_corrected,digits=2)))
      Table2[term,paste0(organ,"__",resp)]<-paste0(teststat,";",p)
    }
  }
}
write.table(Table2,file="tables/Table2_wUF_stats.txt",sep='\t',row.names=FALSE,col.names=TRUE)

####### Combine all LS means into a single dataframe #####
LSM.all<-rbind(mutate(LSM.sha.leaf3,Response="Shannon",Organ="leaf",Term=rownames(LSM.sha.leaf3)),
               mutate(LSM.chao1.leaf3,Response="Chao1",Organ="leaf",Term=rownames(LSM.chao1.leaf3)),
               mutate(LSM.wUF1.leaf3,Response="wUF1",Organ="leaf",Term=rownames(LSM.wUF1.leaf3)),
               mutate(LSM.wUF2.leaf3,Response="wUF2",Organ="leaf",Term=rownames(LSM.wUF2.leaf3)),
               mutate(LSM.wUF3.leaf3,Response="wUF3",Organ="leaf",Term=rownames(LSM.wUF3.leaf3)),
               mutate(LSM.sha.root3,Response="Shannon",Organ="root",Term=rownames(LSM.sha.root3)),
               mutate(LSM.chao1.root3,Response="Chao1",Organ="root",Term=rownames(LSM.chao1.root3)),
               mutate(LSM.wUF1.root3,Response="wUF1",Organ="root",Term=rownames(LSM.wUF1.root3)),
               mutate(LSM.wUF2.root3,Response="wUF2",Organ="root",Term=rownames(LSM.wUF2.root3)),
               mutate(LSM.wUF3.root3,Response="wUF3",Organ="root",Term=rownames(LSM.wUF3.root3)))
rm(LSM.sha.leaf3,LSM.chao1.leaf3,LSM.wUF1.leaf3,LSM.wUF2.leaf3,LSM.wUF3.leaf3)
rm(LSM.sha.root3,LSM.chao1.root3,LSM.wUF1.root3,LSM.wUF2.root3,LSM.wUF3.root3)
save(LSM.all,file="LMMs/LSmeans.RData")


####### Effect size: change in Chao1 richness between genotypes #######
filter(LSM.all,grepl("Genotype ",Term) & Organ=='leaf' & Response=='Chao1') %>% arrange(Estimate)
# smallest LS mean = 967.8565 for MAH, largest = 1168.6981 for JAM

(1168.6981-967.8565)/967.8565 # increase of 20.8% from MAH to JAM genotypes
(1168.6981-976.7557)/976.7557 # increase of 19.7% from MIL to JAM genotypes

####### Fig. S6: alpha diversity ~ Year #######
# Chao1 diversity
pdf(file="plots/Fig_S7a_LSM_Chao1_by_Year.pdf",height=9)
filter(LSM.all,grepl("^Harvested ", Term), Response=='Chao1') %>%
  ggplot(.,aes(x=Harvested, y=Estimate,color=Harvested))+
  facet_wrap(~Organ,ncol=2,scales="free")+
  geom_errorbar(size=2,aes(ymin=Estimate-SE, ymax=Estimate+SE), width=.1,position=position_dodge(width=0.3))+
  geom_point(size=4,position=position_dodge(width=0.3))+
  ylab("Chao1 estimated richness\n")+xlab("\nYear harvested")+
  scale_colour_manual(values=c("red","blue"),guide=FALSE)+
  theme_classic()+
  theme(axis.title.x=element_text(size=28,face="bold"),axis.text.x=element_text(size=24,face="bold"))+
  theme(axis.title.y=element_text(size=28,face="bold"),axis.text.y=element_text(size=24,face="bold"))+
  theme(strip.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))
dev.off()

# Shannon diversity
pdf(file="plots/Fig_S7b_LSM_Shannon_by_Year.pdf",height=9)
filter(LSM.all,grepl("^Harvested ", Term), Response=='Shannon') %>%
  ggplot(.,aes(x=Harvested, y=Estimate,color=Harvested))+
  facet_wrap(~Organ,ncol=2,scales="free")+
  geom_errorbar(size=2,aes(ymin=Estimate-SE, ymax=Estimate+SE), width=.1,position=position_dodge(width=0.3))+
  geom_point(size=4,position=position_dodge(width=0.3))+
  ylab("Shannon effective diversity\n")+xlab("\nYear harvested")+
  scale_colour_manual(values=c("red","blue"),guide=FALSE)+
  scale_y_continuous(labels=function(y) format(exp(y),digits=3),breaks=c(log(35),log(45),log(55),log(430),log(460),log(490)))+ # convert from Shannon entropy to effective Shannon diversity
  theme_classic()+
  theme(axis.title.x=element_text(size=28,face="bold"),axis.text.x=element_text(size=24,face="bold"))+
  theme(axis.title.y=element_text(size=28,face="bold"),axis.text.y=element_text(size=24,face="bold"))+
  theme(strip.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))
dev.off()

####### Fig. 2a: Boxplots: Alpha diversity by Site #######
pdf(file="plots/Fig_2a_Boxplot_Chao1_by_Site.pdf",width=9,height=9)
update_geom_defaults("point", list(colour = NULL))
rbind(mutate(leaf3.smd,Organ='leaf'),
      mutate(root3.smd,Organ='root')) %>%
  ggplot(.,aes(x=Organ,y=Chao1,colour=Site))+
  geom_boxplot(size=1.5,outlier.colour=NULL,position=position_dodge(w=0.9),outlier.size=3)+
  ylab("Chao1 estimated richness\n")+
  scale_color_manual(values=sitePalette)+
  theme_classic()+
  theme(axis.text.x = element_text(size=32,face="bold"))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=32,face="bold"))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=28,face="bold"))+
  theme(axis.text.y=element_text(size=28,face="bold"))+
  theme(legend.background = element_rect(fill="gray95", size=.5))+
  theme(legend.key.height=unit(2,"line"),legend.key.width=unit(2,"line"))
update_geom_defaults("point", list(colour = "black"))
dev.off()

pdf(file="plots/Fig_2b_Boxplot_Shannon_by_Site.pdf",width=9,height=9)
update_geom_defaults("point", list(colour = NULL))
rbind(mutate(leaf3.smd,Organ='leaf'),
      mutate(root3.smd,Organ='root')) %>%
  ggplot(.,aes(x=Organ,y=expShannon,colour=Site))+
  geom_boxplot(size=1.5,outlier.colour=NULL,position=position_dodge(w=0.9),outlier.size=3)+
  ylab("Effective Shannon diversity")+
  scale_color_manual(values=sitePalette)+
  theme_classic()+
  theme(axis.text.x = element_text(size=32,face="bold"))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=32,face="bold"))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=28,face="bold"))+
  theme(axis.text.y=element_text(size=28,face="bold"))+
  theme(legend.background = element_rect(fill="gray95", size=.5))+
  theme(legend.key.height=unit(2,"line"),legend.key.width=unit(2,"line"))
update_geom_defaults("point", list(colour = "black"))
dev.off()

####### Figure 2c: Chao1 Year x Site #######
pdf(file="plots/Fig_2c_LSM_Chao1_YxS_leaf.pdf",width=9,height=6)
filter(LSM.all, grepl("Site:Harvested",Term), Response=='Chao1', Organ=='leaf') %>%
  ggplot(.,aes(x=Site,y=Estimate,color=Harvested))+
  geom_point(size=4,position=position_dodge(w=0.3))+
  geom_errorbar(aes(ymin=Estimate-SE,ymax=Estimate+SE),width=0,size=2,position=position_dodge(w=0.3))+
  scale_color_manual(values=c('red','blue'),name='Year')+
  theme_classic()+
  ggtitle("Leaves")+
  ylab("Chao1 richness")+
  theme(plot.title=element_text(size=36,face="bold",colour="forest green"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,vjust=0.8,hjust=0.9,face="bold",colour=sitePalette))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=28,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

pdf(file="plots/Fig_2c_LSM_Chao1_YxS_root.pdf",width=9,height=6)
filter(LSM.all, grepl("Site:Harvested",Term), Response=='Chao1', Organ=='root') %>%
  ggplot(.,aes(x=Site,y=Estimate,color=Harvested))+
  geom_point(size=4,position=position_dodge(w=0.3))+
  geom_errorbar(aes(ymin=Estimate-SE,ymax=Estimate+SE),width=0,size=2,position=position_dodge(w=0.3))+
  scale_color_manual(values=c('red','blue'),name='Year')+
  theme_classic()+
  ggtitle("Roots")+
  ylab("Chao1 richness")+
  theme(plot.title=element_text(size=36,face="bold",colour="grey"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,vjust=0.8,hjust=0.9,face="bold",colour=sitePalette))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=28,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

####### Fig. 2d: Leaf&Root: weighted UniFrac LSM ordination ~ Year #######
pdf(file="plots/Fig_2d_LSMordination_byYear_leaf.pdf",width=9,height=9)
filter(LSM.all, grepl("Harvested ",Term), Response%in%c('wUF1','wUF2'), Organ=='leaf') %>%
  reshape(v.names=c('Estimate','SE'),timevar='Response',direction='wide',idvar=c('Harvested')) %>%
  plyr::rename(replace=c('Estimate.wUF1'='wUF1','Estimate.wUF2'='wUF2')) %>%
  ggplot(.,aes(x=wUF1,y=wUF2,color=Harvested))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=wUF2-SE.wUF2,ymax=wUF2+SE.wUF2),width=0.00,size=2)+
  geom_errorbarh(aes(xmin=wUF1-SE.wUF1,xmax=wUF1+SE.wUF1),height=0.00,size=2)+
  scale_color_manual(values=c('red','blue'),name='Year')+
  ylab("PCo2 [11.2%]")+xlab("PCo1 [45.3%]")+
  theme_classic()+
  theme(plot.title = element_text(size=40, face="bold"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=32,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=32,face="bold"))+
  theme(legend.title= element_text(size=34),legend.text=element_text(size=30,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(1.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

pdf(file="plots/Fig_2d_LSMordination_byYear_root.pdf",width=9,height=9)
filter(LSM.all, grepl("Harvested ",Term), Response%in%c('wUF1','wUF2'), Organ=='root') %>%
  reshape(v.names=c('Estimate','SE'),timevar='Response',direction='wide',idvar=c('Harvested')) %>%
  plyr::rename(replace=c('Estimate.wUF1'='wUF1','Estimate.wUF2'='wUF2')) %>%
  ggplot(.,aes(x=wUF1,y=wUF2,color=Harvested))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=wUF2-SE.wUF2,ymax=wUF2+SE.wUF2),width=0.00,size=2)+
  geom_errorbarh(aes(xmin=wUF1-SE.wUF1,xmax=wUF1+SE.wUF1),height=0.00,size=2)+
  scale_color_manual(values=c('red','blue'),name='Year')+
  ylab("PCo2 [16.1%]")+xlab("PCo1 [29.6%]")+
  ggtitle("Roots\n")+theme_classic()+
  theme(plot.title = element_text(size=40, face="bold"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=32,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=32,face="bold"))+
  theme(legend.title= element_text(size=34),legend.text=element_text(size=30,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(1.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

####### Fig. 5c: Leaf&Root: weighted UniFrac LSM ordination ~ Genotype #######
pdf(file="plots/Fig5c_LSMordination_byGenotype_leaf.pdf",width=9,height=9)
filter(LSM.all, grepl("Genotype ",Term), Response%in%c('wUF1','wUF3'), Organ=='leaf') %>%
  reshape(v.names=c('Estimate','SE'),timevar='Response',direction='wide',idvar=c('Genotype')) %>%
  plyr::rename(replace=c('Estimate.wUF1'='wUF1','Estimate.wUF3'='wUF3')) %>%
  ggplot(.,aes(x=wUF1,y=wUF3,color=Genotype))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=wUF3-SE.wUF3,ymax=wUF3+SE.wUF3),width=0.00,size=2)+
  geom_errorbarh(aes(xmin=wUF1-SE.wUF1,xmax=wUF1+SE.wUF1),height=0.00,size=2)+
  scale_color_manual(values=popPalette)+
  ylab("PCo3 [6.2%]")+xlab("PCo1 [45.2%]")+
  theme_classic()+
  theme(plot.title = element_text(size=40, face="bold"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=32,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=32,face="bold"))+
  theme(legend.title= element_text(size=34),legend.text=element_text(size=30,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(1.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

pdf(file="plots/Fig5c_LSMordination_byGenotype_root.pdf",width=9,height=9)
filter(LSM.all, grepl("Genotype ",Term), Response%in%c('wUF1','wUF2'), Organ=='root') %>%
  reshape(v.names=c('Estimate','SE'),timevar='Response',direction='wide',idvar=c('Genotype')) %>%
  plyr::rename(replace=c('Estimate.wUF1'='wUF1','Estimate.wUF2'='wUF2')) %>%
  ggplot(.,aes(x=wUF1,y=wUF2,color=Genotype))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=wUF2-SE.wUF2,ymax=wUF2+SE.wUF2),width=0.00,size=2)+
  geom_errorbarh(aes(xmin=wUF1-SE.wUF1,xmax=wUF1+SE.wUF1),height=0.00,size=2)+
  scale_color_manual(values=popPalette)+
  scale_x_continuous(breaks=c(-0.0125,-0.0075,-0.0025))+
  ylab("PCo2 [16.1%]")+xlab("PCo1 [29.6%]")+
  ggtitle("Roots\n")+theme_classic()+
  theme(plot.title = element_text(size=40, face="bold"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=32,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=32,face="bold"))+
  theme(legend.title= element_text(size=34),legend.text=element_text(size=30,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(1.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

####### Fig. 6b: Leaf- weighted UniFrac LSM ordination ~ Genotype x Site #######

pdf(file="plots/Fig6b_wUF1_LSM_GxS_leaf.pdf",width=9,height=6)
filter(LSM.all, grepl("Genotype:Site",Term), Response=='wUF1', Organ=='leaf') %>%
ggplot(.,aes(x=Site,y=Estimate,color=Genotype))+
  geom_point(size=4,position=position_dodge(w=0.45))+
  geom_errorbar(aes(ymin=Estimate-SE,ymax=Estimate+SE),width=0,size=2,position=position_dodge(w=0.45))+
  scale_color_manual(values=popPalette)+
  theme_classic()+
  ylab("PCo1 [46.1%]")+
  coord_cartesian(ylim=c(-0.05,0.10))+ # to get the whole range to show up on vertical axis
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,vjust=0.8,hjust=0.9,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=28,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()
  
####### Figure: Shannon diversity ~ GxS, leaf ####### 
pdf(file="plots/Fig_6c_LSM_Shannon_GxS_leaf.pdf",width=9,height=6)
filter(LSM.all, grepl("Genotype:Site",Term), Response=='Shannon', Organ=='leaf') %>%
  ggplot(.,aes(x=Site,y=Estimate,color=Genotype))+
  geom_errorbar(size=2,aes(ymin=Estimate-SE, ymax=Estimate+SE), width=.1,position=position_dodge(width=0.45))+
  geom_point(size=4,position=position_dodge(width=0.45))+
  ylab("Shannon diversity\n")+xlab("\nSite")+
  scale_color_manual(values=popPalette)+
  theme_classic()+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,vjust=0.8,hjust=0.9,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=28,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
 dev.off()

####### Figure 4c: Chao1 ~ Age x Site #######

pdf(file="plots/Fig_4c_LSM_Chao1_AxS_leaf.pdf",width=9,height=6)
filter(LSM.all, grepl("Site:Age",Term), Response=='Chao1', Organ=='leaf') %>%
  ggplot(.,aes(x=Site,y=Estimate,color=Age))+
  geom_point(size=4,position=position_dodge(w=0.3))+
  geom_errorbar(aes(ymin=Estimate-SE,ymax=Estimate+SE),width=0,size=2,position=position_dodge(w=0.3))+
  scale_color_manual(values=agePalette)+
  theme_classic()+
  ggtitle("Leaves")+
  ylab("Chao1 richness")+
  theme(plot.title=element_text(size=36,face="bold",colour="forest green"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,vjust=0.8,hjust=0.9,face="bold",colour=sitePalette))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=28,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

pdf(file="plots/Fig_4c_LSM_Chao1_AxS_root.pdf",width=9,height=6)
filter(LSM.all, grepl("Site:Age",Term), Response=='Chao1', Organ=='root') %>%
  ggplot(.,aes(x=Site,y=Estimate,color=Age))+
  geom_point(size=4,position=position_dodge(w=0.3))+
  geom_errorbar(aes(ymin=Estimate-SE,ymax=Estimate+SE),width=0,size=2,position=position_dodge(w=0.3))+
  scale_color_manual(values=agePalette)+
  theme_classic()+
  ggtitle("Roots")+
  ylab("Chao1 richness")+
  theme(plot.title=element_text(size=36,face="bold",colour="grey"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,vjust=0.8,hjust=0.9,face="bold",colour=sitePalette))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=28,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

####### Fig. 5a: LSM Chao1 ~ Genotype #######
pdf(file="plots/Fig_5a_LSM_Chao1_byGenotype.pdf",height=9,width=12)
filter(LSM.all,grepl("^Genotype ", Term), Response=='Chao1') %>%
  ggplot(.,aes(x=Genotype,y=Estimate,color=Genotype))+
  facet_wrap(~Organ,scales="free")+
  geom_errorbar(size=2,aes(ymin=Estimate-SE, ymax=Estimate+SE), width=.1,position=position_dodge(width=0.3))+
  geom_point(size=4,position=position_dodge(width=0.3))+
  ylab("Chao1 estimated richness\n")+xlab("\nGenotype")+
  scale_colour_manual(values=popPalette,guide=FALSE)+
  theme_classic()+theme(panel.margin=unit(2,"lines"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=28,face="bold",angle=35,vjust=0.6))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=28,face="bold"))+
  theme(strip.background=element_rect(fill="gray90",color="gray90"),strip.text=element_text(size=36,face="bold"))
dev.off()

####### Fig. S12: LSM Shannon ~ Genotype #######
pdf(file="plots/Fig_S12_LSM_Shannon_byGenotype.pdf",height=5,width=11)
filter(LSM.all,grepl("^Genotype ", Term), Response=='Shannon') %>%
  ggplot(.,aes(x=Genotype,y=exp(Estimate),color=Genotype))+
  facet_wrap(~Organ,ncol=2,scales="free")+
  geom_errorbar(size=3, aes(ymin=exp(Estimate-SE), ymax=exp(Estimate+SE)), width=.1,position=position_dodge(width=0.3))+
  geom_point(size=6,position=position_dodge(width=0.3))+
  ylab("Shannon\neffective diversity\n")+xlab("\nGenotype")+
  scale_colour_manual(values=popPalette,guide=FALSE)+
  theme_classic()+
  theme(strip.text.x=element_text(size=28,face="bold"),strip.background=element_rect(color="gray90",fill="gray90"))+
  theme(axis.title.x=element_text(size=28,face="bold"),axis.text.x=element_text(size=23,face="bold"))+
  theme(axis.title.y=element_text(size=26,face="bold"),axis.text.y=element_text(size=25,face="bold"))
dev.off()

####### Fig. 4a: LSM Chao1 ~ Age #######
# Chao1 diversity: both tissues
pdf(file="plots/Fig_4a_LSM_Chao1_byAge.pdf",height=9,width=15)
filter(LSM.all,grepl("^Age ", Term), Response=='Chao1') %>%
  ggplot(.,aes(x=Age,y=Estimate,color=Age))+
  facet_wrap(~Organ,scales="free")+
  geom_errorbar(size=3,aes(ymin=Estimate-SE, ymax=Estimate+SE), width=.05,position=position_dodge(width=0.3))+
  geom_point(size=7,position=position_dodge(width=0.3))+
  ylab("Chao1 estimated richness\n")+xlab("Plant age (years)")+
  scale_colour_manual(values=agePalette,guide=FALSE)+
  theme(panel.margin = unit(1.5, "lines"))+
  theme_classic()+guides(color=FALSE)+
  theme(panel.margin=unit(2,"lines"))+
  theme(axis.title.x=element_text(size=38,face="bold"),axis.text.x=element_text(size=36,face="bold"))+
  theme(axis.title.y=element_text(size=38,face="bold",lineheight=0.8),axis.text.y=element_text(size=34,face="bold"))+
  theme(strip.text.x=element_text(size=38,face="bold"),strip.background=element_rect(fill="gray90",color="gray90"))
dev.off()

####### Save image #######
save.image(paste0("LMMs/image_",date(),".RData"))
