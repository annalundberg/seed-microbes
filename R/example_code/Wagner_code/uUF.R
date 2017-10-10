### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
##### Parallel analysis: Linear mixed models of beta diversity using the unweighted UniFrac metric
####### Clear workspace ########
rm(list=ls())

####### Load source file #######
source('ecotypes_source.R')

####### Session Info #######
sessionInfo()
"R version 3.2.3 (2015-12-10)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tidyr_0.4.1                VennDiagram_1.6.17         futile.logger_1.4.1        BiocParallel_1.4.3        
[5] reshape2_1.4.1             lmerTest_2.0-30            lme4_1.1-11                Matrix_1.2-4              
[9] scales_0.4.0               ggplot2_2.1.0              vegan_2.3-4                lattice_0.20-33           
[13] permute_0.9-0              doParallel_1.0.10          iterators_1.0.8            foreach_1.4.3             
[17] dplyr_0.4.3                plyr_1.8.3                 genefilter_1.52.1          Biostrings_2.38.4         
[21] XVector_0.10.0             ape_3.4                    mapdata_2.2-6              maps_3.1.0                
[25] DESeq2_1.10.1              RcppArmadillo_0.6.600.4.0  Rcpp_0.12.3                SummarizedExperiment_1.0.2
[29] Biobase_2.30.0             GenomicRanges_1.22.4       GenomeInfoDb_1.6.3         IRanges_2.4.8             
[33] S4Vectors_0.8.11           BiocGenerics_0.16.1        phyloseq_1.14.0           

loaded via a namespace (and not attached):
[1] splines_3.2.3        Formula_1.2-1        assertthat_0.1       latticeExtra_0.6-28  RSQLite_1.0.0       
[6] digest_0.6.9         chron_2.3-47         RColorBrewer_1.1-2   minqa_1.2.4          colorspace_1.2-6    
[11] XML_3.98-1.4         zlibbioc_1.16.0      xtable_1.8-2         annotate_1.48.0      mgcv_1.8-12         
[16] lazyeval_0.1.10      nnet_7.3-12          survival_2.38-3      RJSONIO_1.3-0        magrittr_1.5        
[21] nlme_3.1-126         MASS_7.3-45          foreign_0.8-66       tools_3.2.3          data.table_1.9.6    
[26] stringr_1.0.0        munsell_0.4.3        locfit_1.5-9.1       cluster_2.0.3        AnnotationDbi_1.32.3
[31] lambda.r_1.1.7       compiler_3.2.3       ade4_1.7-4           nloptr_1.0.4         biom_0.3.12         
[36] igraph_1.0.1         labeling_0.3         gtable_0.2.0         codetools_0.2-14     multtest_2.26.0     
[41] DBI_0.3.1            R6_2.1.2             gridExtra_2.2.1      Hmisc_3.17-2         futile.options_1.0.0
[46] stringi_1.0-1        geneplotter_1.48.0   rpart_4.1-10         acepack_1.3-3.3     "
####### Load data: variance-stabilizing transformed Phyloseq objects for roots and leaves at 3 sites ######

load("intermediate_data/phylo_leaf3_withEndog_vst.RData")
load("intermediate_data/phylo_root3_withEndog_vst.RData")
####### Register parallel backend #######
registerDoParallel(cores=4)
####### UNweighted UniFrac and PCoA: separately for leaf and root datasets #######
# replace negative values with 0s just for distance calculations
uUF.root3.withEndog.vst<-UniFrac(
  transform_sample_counts(root3.withEndog.vst,function(x) x<-ifelse(x<0,0,x)),
  weighted=FALSE,parallel=TRUE)
uUF.leaf3.withEndog.vst<-UniFrac(
  transform_sample_counts(leaf3.withEndog.vst,function(x) x<-ifelse(x<0,0,x)),
  weighted=FALSE,parallel=TRUE)

save(uUF.leaf3.withEndog.vst,file="uUF/uUF_leaf3_wEndog_vst.RData")
save(uUF.root3.withEndog.vst,file="uUF/uUF_root3_wEndog_vst.RData")

cap.uUF.root3.withEndog.vst<-capscale(uUF.root3.withEndog.vst~1,data=as(sample_data(root3.withEndog.vst),'data.frame'))
cap.uUF.leaf3.withEndog.vst<-capscale(uUF.leaf3.withEndog.vst~1,data=as(sample_data(leaf3.withEndog.vst),'data.frame'))

## Get inertia for top 3 PCoA axes: leaf ##
cap.uUF.leaf3.withEndog.vst$CA$eig[1:3]/sum(cap.uUF.leaf3.withEndog.vst$CA$eig) 
"      MDS1       MDS2       MDS3 
0.18686611 0.12875009 0.03706156 "
cap.uUF.root3.withEndog.vst$CA$eig[1:3]/sum(cap.uUF.root3.withEndog.vst$CA$eig) 
"      MDS1       MDS2       MDS3 
0.30832924 0.09191178 0.05417460 "

####### Fig. S5 Scree plots: unweighted UniFrac #######
pdf(file="plots/Fig_S5c_Supp_screeplot_uUF_leaf.pdf")
barplot(cap.uUF.leaf3.withEndog.vst$CA$eig[1:20]/sum(cap.uUF.leaf3.withEndog.vst$CA$eig),
        main="Leaves: unweighted UniFrac",xlab="Principal coordinate",ylab="Proportion Variance",
        cex.lab=2,cex.main=2.5,cex.axis=1.5,font.axis=2,font.main=2,font.lab=2,xaxt='n')
dev.off()
pdf(file="plots/Fig_S5d_screeplot_uUF_root.pdf")
barplot(cap.uUF.root3.withEndog.vst$CA$eig[1:20]/sum(cap.uUF.root3.withEndog.vst$CA$eig),
        main="Roots: unweighted UniFrac",xlab="Principal coordinate",ylab="Proportion Variance",
        cex.lab=2,cex.main=2.5,cex.axis=1.5,font.axis=2,font.main=2,font.lab=2,xaxt='n')
dev.off()
####### How much variation is explained by the top 3 PCo axes? #######
sink("ordination_top3_cumulative_PVEs.txt")
print("cumulative percent variance explained by top 3 unweighted UniFrac PCo:")
print("unweighted UniFrac, roots:")
sum(cap.uUF.root3.withEndog.vst$CA$eig[1:3])/sum(cap.uUF.root3.withEndog.vst$CA$eig)  # 0.4544
print("unweighted UniFrac, leaves:")
sum(cap.uUF.leaf3.withEndog.vst$CA$eig[1:3])/sum(cap.uUF.leaf3.withEndog.vst$CA$eig)  # 0.3527
print("Individual percent variance explained by top 3 unweighted UniFrac PCo:")
print("unweighted UniFrac, roots:")
(cap.uUF.root3.withEndog.vst$CA$eig[1:3])/sum(cap.uUF.root3.withEndog.vst$CA$eig)  
print("unweighted UniFrac, leaves:")
(cap.uUF.leaf3.withEndog.vst$CA$eig[1:3])/sum(cap.uUF.leaf3.withEndog.vst$CA$eig) 
sink()

####### Fig. unweighted UniFrac Ordination ~ Site #######
pdf("plots/Fig_S4a_Ordination_uUF1_2_Site_leaf.pdf",width=9,height=9)
plot_ordination(leaf3.withEndog.vst,cap.uUF.leaf3.withEndog.vst,type="samples",axes=1:2,color="Site")+
  scale_colour_manual(values=sitePalette)+
  geom_point(size=4,alpha=1)+
  xlab("unweighted UniFrac\nPCo1 [18.7%]")+ylab("unweighted UniFrac\nPCo2 [12.9%]")+
  ggtitle("Leaves")+theme_classic()+
  theme(plot.title = element_text(size=44, face="bold",color="forest green"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=40),legend.text=element_text(size=36,face="bold"))+
  theme(legend.key.height=unit(2.5,"lines"),legend.key.width=unit(2,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

pdf("plots/Fig_S4b_Ordination_uUF1_2_Site_root.pdf",width=9,height=9)
plot_ordination(root3.withEndog.vst,cap.uUF.root3.withEndog.vst,type="samples",axes=1:2,color="Site")+
  scale_colour_manual(values=sitePalette,guide=FALSE)+
  geom_point(size=4,alpha=1)+
  xlab("unweighted UniFrac\nPCo1 [30.8%]")+ylab("unweighted UniFrac\nPCo2 [9.2%]")+
  ggtitle("Roots")+theme_classic()+
  theme(plot.title = element_text(size=44, face="bold",color="grey"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))
dev.off()

####### Save major PCoA axes and Alpha diversity metrics for use in LMMs: #######
leaf3.smd.uUF.withEndog<-as(sample_data(leaf3.withEndog.vst),'data.frame') %>%
  mutate(SampleID=row.names(.)) %>%
  merge(.,as.data.frame(cap.uUF.leaf3.withEndog.vst$CA$u[,1:3]),by.x='SampleID',by.y="row.names") %>%
  plyr::rename(replace=c('MDS1'='uUF1','MDS2'='uUF2','MDS3'='uUF3'))

root3.smd.uUF.withEndog<-as(sample_data(root3.withEndog.vst),'data.frame') %>%
  mutate(SampleID=row.names(.)) %>%
  merge(.,as.data.frame(cap.uUF.root3.withEndog.vst$CA$u[,1:3]),by.x='SampleID',by.y="row.names") %>%
  plyr::rename(replace=c('MDS1'='uUF1','MDS2'='uUF2','MDS3'='uUF3'))

save(leaf3.smd.uUF.withEndog,file="uUF/smd_leaf3_withEndog.RData")
save(root3.smd.uUF.withEndog,file="uUF/smd_root3_withEndog.RData")

####### Remove endogenous plants for model fitting #######
leaf3.smd.uUF<-filter(leaf3.smd.uUF.withEndog,Age!='endog')
root3.smd.uUF<-filter(root3.smd.uUF.withEndog,Age!='endog')

####### Check residuals and R^2 for unweighted UniFrac models and save LS means for plotting #######
# unweighted UniFrac PCo1: leaves
leaf3.uUF1<-lmer(uUF1~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=leaf3.smd.uUF,REML=TRUE)
plot(residuals(leaf3.uUF1)~fitted(leaf3.uUF1))
qqnorm(residuals(leaf3.uUF1)); qqline(residuals(leaf3.uUF1))  
LSM.uUF1.leaf3<-as.data.frame(lmerTest::lsmeans(leaf3.uUF1)[1])
colnames(LSM.uUF1.leaf3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(leaf3.uUF1) # 0.85085
rm(leaf3.uUF1)

# unweighted UniFrac PCo2: leaves
leaf3.uUF2<-lmer(uUF2~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=leaf3.smd.uUF,REML=TRUE)
plot(residuals(leaf3.uUF2)~fitted(leaf3.uUF2))
qqnorm(residuals(leaf3.uUF2)); qqline(residuals(leaf3.uUF2)) 
LSM.uUF2.leaf3<-as.data.frame(lmerTest::lsmeans(leaf3.uUF2)[1])
colnames(LSM.uUF2.leaf3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(leaf3.uUF2) # 0.80380
rm(leaf3.uUF2)

# unweighted UniFrac PCo3: leaves
leaf3.uUF3<-lmer(uUF3~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=leaf3.smd.uUF,REML=TRUE)
plot(residuals(leaf3.uUF3)~fitted(leaf3.uUF3))
qqnorm(residuals(leaf3.uUF3)); qqline(residuals(leaf3.uUF3))
LSM.uUF3.leaf3<-as.data.frame(lmerTest::lsmeans(leaf3.uUF3)[1])
colnames(LSM.uUF3.leaf3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(leaf3.uUF3) # 0.8697594
rm(leaf3.uUF3)

# unweighted UniFrac PCo1: roots
root3.uUF1<-lmer(uUF1~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=root3.smd.uUF,REML=TRUE)
plot(residuals(root3.uUF1)~fitted(root3.uUF1)) 
qqnorm(residuals(root3.uUF1)); qqline(residuals(root3.uUF1)) 
LSM.uUF1.root3<-as.data.frame(lmerTest::lsmeans(root3.uUF1)[1])
colnames(LSM.uUF1.root3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(root3.uUF1) # 0.9784
rm(root3.uUF1)

# unweighted UniFrac PCo2: roots
root3.uUF2<-lmer(uUF2~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=root3.smd.uUF,REML=TRUE)
plot(residuals(root3.uUF2)~fitted(root3.uUF2)) 
qqnorm(residuals(root3.uUF2)); qqline(residuals(root3.uUF2)) 
LSM.uUF2.root3<-as.data.frame(lmerTest::lsmeans(root3.uUF2)[1])
colnames(LSM.uUF2.root3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(root3.uUF2) # 0.9575643
rm(root3.uUF2)

# unweighted UniFrac PCo3: roots
root3.uUF3<-lmer(uUF3~Genotype*Site + Age*Site + Harvested*Site + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=root3.smd.uUF,REML=TRUE)
plot(residuals(root3.uUF3)~fitted(root3.uUF3))
qqnorm(residuals(root3.uUF3)); qqline(residuals(root3.uUF3)) 
LSM.uUF3.root3<-as.data.frame(lmerTest::lsmeans(root3.uUF3)[1])
colnames(LSM.uUF3.root3)<-c("Genotype","Site","Age","Harvested","Estimate","SE","DF","t","lowerCI","upperCI","P_uncorrected")
r2.LMM(root3.uUF3) # 0.91394
rm(root3.uUF3)

####### Table S5: unweighted UniFrac results #######
# use "LMMstats" function to control for MiSeq run and sequencing depth (logObs)
stats.uUF.leaf3<-LMMstats(leaf3.smd.uUF,resp="uUF1",reml=TRUE) %>%
  rbind(.,LMMstats(leaf3.smd.uUF,resp="uUF2",reml=TRUE)) %>%
  rbind(.,LMMstats(leaf3.smd.uUF,resp="uUF3",reml=TRUE)) %>%
  Pcorrect(method='holm')

stats.uUF.root3<-LMMstats(root3.smd.uUF,resp="uUF1",reml=TRUE) %>%
  rbind(.,LMMstats(root3.smd.uUF,resp="uUF2",reml=TRUE)) %>%
  rbind(.,LMMstats(root3.smd.uUF,resp="uUF3",reml=TRUE)) %>%
  Pcorrect(method='holm')

TableUUF<-data.frame(row.names=levels(stats.uUF.leaf3$Term))
TableUUF$Term<-rownames(TableUUF)
for (h in 1:2){
  organ<-c('leaf','root')[h]
  stats<-get(c("stats.uUF.leaf3","stats.uUF.root3")[h])
  for (i in 1:length(levels(stats$Response))){
    resp<-levels(stats$Response)[i]
    for (j in 1:length(levels(stats$Term))){
      term<-levels(stats$Term)[j]
      subTableUUF<-subset(stats,Response==resp&Term==term)
      teststat<-ifelse(term%in%c('Block','Line','Plate'),
                       paste0("ChiSq",subTableUUF$df,"=",format(subTableUUF$F.or.ChiSq,digits=3)),
                       paste0("F",subTableUUF$df,",",ceiling(subTableUUF$ddf),"=",format(subTableUUF$F.or.ChiSq,digits=3)))
      p<-ifelse(subTableUUF$P_corrected==0,"P<3e-16",paste0("P=",format(subTableUUF$P_corrected,digits=2)))
      TableUUF[term,paste0(organ,"__",resp)]<-paste0(teststat,";",p)
    }
  }
}
write.table(TableUUF,file="tables/Table_S5_UUFstats.txt",sep='\t',row.names=FALSE,col.names=TRUE)

####### Combine all LS means into a single dataframe #####
LSM.uUF.all<-rbind(mutate(LSM.uUF1.leaf3,Response="uUF1",Organ="leaf",Term=rownames(LSM.uUF1.leaf3)),
               mutate(LSM.uUF2.leaf3,Response="uUF2",Organ="leaf",Term=rownames(LSM.uUF2.leaf3)),
               mutate(LSM.uUF3.leaf3,Response="uUF3",Organ="leaf",Term=rownames(LSM.uUF3.leaf3)),
               mutate(LSM.uUF1.root3,Response="uUF1",Organ="root",Term=rownames(LSM.uUF1.root3)),
               mutate(LSM.uUF2.root3,Response="uUF2",Organ="root",Term=rownames(LSM.uUF2.root3)),
               mutate(LSM.uUF3.root3,Response="uUF3",Organ="root",Term=rownames(LSM.uUF3.root3)))
rm(LSM.uUF1.leaf3,LSM.uUF2.leaf3,LSM.uUF3.leaf3)
rm(LSM.uUF1.root3,LSM.uUF2.root3,LSM.uUF3.root3)
save(LSM.uUF.all,file="uUF/uUF_LSmeans.RData")

####### Figure S13 : Leaf&Root: unweighted UniFrac LSM ordination ~ Genotype #######
pdf(file="plots/Fig_S13a_uUF_LSMordination_byGenotype_leaf.pdf",width=9,height=9)
filter(LSM.uUF.all, grepl("Genotype ",Term), Response%in%c('uUF1','uUF3'), Organ=='leaf') %>%
  reshape(v.names=c('Estimate','SE'),timevar='Response',direction='wide',idvar=c('Genotype')) %>%
  plyr::rename(replace=c('Estimate.uUF1'='uUF1','Estimate.uUF3'='uUF3')) %>%
  ggplot(.,aes(x=uUF1,y=uUF3,color=Genotype))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=uUF3-SE.uUF3,ymax=uUF3+SE.uUF3),width=0.00,size=2)+
  geom_errorbarh(aes(xmin=uUF1-SE.uUF1,xmax=uUF1+SE.uUF1),height=0.00,size=2)+
  scale_color_manual(values=popPalette)+
  ylab("unweighted UniFrac\nPCo3 [3.7%]")+xlab("unweighted UniFrac\nPCo1 [18.7%]")+
  ggtitle("Leaves\n")+theme_classic()+
  theme(plot.title = element_text(size=40, face="bold",color="forest green"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=25,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=25,face="bold"))+
  theme(legend.title= element_text(size=34),legend.text=element_text(size=30,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(1.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

pdf(file="plots/Fig_S13b_uUF_LSMordination_byGenotype_root.pdf",width=9,height=9)
filter(LSM.uUF.all, grepl("Genotype ",Term), Response%in%c('uUF1','uUF2'), Organ=='root') %>%
  reshape(v.names=c('Estimate','SE'),timevar='Response',direction='wide',idvar=c('Genotype')) %>%
  plyr::rename(replace=c('Estimate.uUF1'='uUF1','Estimate.uUF2'='uUF2')) %>%
  ggplot(.,aes(x=uUF1,y=uUF2,color=Genotype))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=uUF2-SE.uUF2,ymax=uUF2+SE.uUF2),width=0.00,size=2)+
  geom_errorbarh(aes(xmin=uUF1-SE.uUF1,xmax=uUF1+SE.uUF1),height=0.00,size=2)+
  scale_color_manual(values=popPalette)+
  scale_x_continuous(breaks=c(-0.0125,-0.0075,-0.0025))+
  ylab("unweighted UniFrac\nPCo2 [9.2%]")+xlab("unweighted UniFrac\nPCo1 [30.8%]")+
  ggtitle("Roots\n")+theme_classic()+
  theme(plot.title = element_text(size=40, face="bold",color="dark grey"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=25,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=25,face="bold"))+
  theme(legend.title= element_text(size=34),legend.text=element_text(size=30,face="bold"))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(1.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

####### Save image #######
save.image(paste0("uUF/image_",date(),".RData"))
