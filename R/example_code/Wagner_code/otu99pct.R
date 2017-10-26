### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
#####  Analysis of OTUs grouped at 99% sequence identity
####### Clear workspace ########
rm(list=ls())

####### Load source file #######
source('ecotypes_source.R')

####### Session info #######
sessionInfo()
"R version 3.2.3 (2015-12-10)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets 
[9] methods   base     

other attached packages:
[1] tidyr_0.4.1                VennDiagram_1.6.17         futile.logger_1.4.1       
[4] BiocParallel_1.4.3         reshape2_1.4.1             lmerTest_2.0-30           
[7] lme4_1.1-11                Matrix_1.2-4               scales_0.4.0              
[10] ggplot2_2.1.0              vegan_2.3-4                lattice_0.20-33           
[13] permute_0.9-0              doParallel_1.0.10          iterators_1.0.8           
[16] foreach_1.4.3              dplyr_0.4.3                plyr_1.8.3                
[19] genefilter_1.52.1          Biostrings_2.38.4          XVector_0.10.0            
[22] ape_3.4                    mapdata_2.2-6              maps_3.1.0                
[25] DESeq2_1.10.1              RcppArmadillo_0.6.600.4.0  Rcpp_0.12.3               
[28] SummarizedExperiment_1.0.2 Biobase_2.30.0             GenomicRanges_1.22.4      
[31] GenomeInfoDb_1.6.3         IRanges_2.4.8              S4Vectors_0.8.11          
[34] BiocGenerics_0.16.1        phyloseq_1.14.0           

loaded via a namespace (and not attached):
[1] splines_3.2.3        Formula_1.2-1        assertthat_0.1       latticeExtra_0.6-28 
[5] RSQLite_1.0.0        digest_0.6.9         chron_2.3-47         RColorBrewer_1.1-2  
[9] minqa_1.2.4          colorspace_1.2-6     XML_3.98-1.4         zlibbioc_1.16.0     
[13] xtable_1.8-2         annotate_1.48.0      mgcv_1.8-12          lazyeval_0.1.10     
[17] nnet_7.3-12          survival_2.38-3      RJSONIO_1.3-0        magrittr_1.5        
[21] nlme_3.1-126         MASS_7.3-45          foreign_0.8-66       tools_3.2.3         
[25] data.table_1.9.6     stringr_1.0.0        munsell_0.4.3        locfit_1.5-9.1      
[29] cluster_2.0.3        AnnotationDbi_1.32.3 lambda.r_1.1.7       compiler_3.2.3      
[33] ade4_1.7-4           nloptr_1.0.4         biom_0.3.12          igraph_1.0.1        
[37] labeling_0.3         gtable_0.2.0         codetools_0.2-14     multtest_2.26.0     
[41] DBI_0.3.1            R6_2.1.2             gridExtra_2.2.1      Hmisc_3.17-2        
[45] futile.options_1.0.0 stringi_1.0-1        geneplotter_1.48.0   rpart_4.1-10        
[49] acepack_1.3-3.3   "

####### Load 99% OTU table with contaminants already removed #######
## Load OTU table into R
otu_table99<-read.table("raw_data/otuTable99.txt",skip=1,sep='\t',header=TRUE,comment.char="")

# Rename first column
colnames(otu_table99)[1]<-'OTU_ID'

# Reformat to Phyloseq-compatible matrix
rownames(otu_table99)<-otu_table99[,1] # store OTU ID numbers as row names
rownames(otu_table99)<-paste("OTU_",rownames(otu_table99),sep="") # add "OTU_" prefix to OTU IDs / row names
otu_table99<-otu_table99[,-1] # remove column of OTU IDs

otu_table99<-otu_table(as.matrix(otu_table99),taxa_are_rows=TRUE) # convert to matrix then convert to Phyloseq "otu table" format

####### Load taxonomic assignments for 99%-binned OTUs #######
tax_table99<-as.matrix(read.table('raw_data/taxAssignments99.txt',sep='\t',header=TRUE))
dimnames(tax_table99)[[1]]<-tax_table99[,1]  # make OTU_IDs the row names
tax_table99<-tax_table99[,-1] # get rid of OTU_ID column
tax_table99<-tax_table(tax_table99) # convert to Phyloseq "tax table" format

####### Make phyloseq object with entire dataset #######
phylo99<-phyloseq(otu_table99,tax_table99)
save(phylo99,file="otu99pct/phyloseq_99pct_nocontam_alldata.RData")
rm(otu_table99,tax_table99)

####### Throw out bad taxa #######
ntaxa(phylo99) # 385275 OTUs
phylo99<-subset_taxa(phylo99, Kingdom!='Unassigned' & Class!='Chloroplast' & Family!='mitochondria')
ntaxa(phylo99) # down to 350845 OTUs

####### Pare down to only samples from the Ecotypes experiment #######
load("intermediate_data/Phyloseq_fieldEco_vst.RData")  # to get the sample metadata
fieldEco_99<-prune_samples(sample_names(fieldEco.vst),phylo99) # reduce 99% dataset to the samples from fullEco
sample_data(fieldEco_99)<-sample_data(fieldEco.vst) # copy sample metadata 
save(fieldEco_99,file="otu99pct/phyloseq_fieldEco_noBadOTUs_noThreshold.RData")

####### Standard thresholding: 25x5 (throw out "non-reproducible" OTUs) #######
threshold1<-kOverA(5,A=25) # set threshold values (require k samples with A reads)
thresh1fieldEco_99<-filter_taxa(fieldEco_99,threshold1,TRUE); ntaxa(thresh1fieldEco_99) # 5862 of 350845 passed the threshold
save(thresh1fieldEco_99, file="otu99pct/phyloseq_fieldEco_noBadOTUs_thresholded.RData")

####### Only keep PAIRED leaf/root samples #######
# Get Plant_IDs that are found in both leaf and root datasets
pairs<-intersect(levels(sample_data(subset_samples(thresh1fieldEco_99,Type=='leaf'))$Plant_ID),
                 levels(sample_data(subset_samples(thresh1fieldEco_99,Type=='root'))$Plant_ID))
length(pairs) # 237 plants with paired data

# Separate into paired-leaf and paired-root datasets and throw out taxa with no observations in each dataset
leaf.paired<-subset_samples(thresh1fieldEco_99,Type=='leaf' & Plant_ID %in% pairs) %>% prune_taxa(taxa_sums(.)>0,.) 
ntaxa(leaf.paired) # 5537 of 5862 OTUs remain
root.paired<-subset_samples(thresh1fieldEco_99,Type=='root' & Plant_ID %in% pairs) %>% prune_taxa(taxa_sums(.)>0,.) 
ntaxa(root.paired) # 5848 of 5862 OTUs remain

####### Get list of OTUs found in both leaves and roots ####### 
shared.otus<-intersect(taxa_names(leaf.paired),taxa_names(root.paired)) # 5523 OTUs are found in both organs
bulksoils<-subset_samples(thresh1fieldEco_99,Type=='soil') %>% prune_taxa(taxa_sums(.)>0,.) # separate out bulk soil data
shared.otus.in.soil<-intersect(taxa_names(bulksoils),shared.otus)
length(shared.otus.in.soil)/length(shared.otus) # 92.6% of OTUs found in both leaves and roots were also found in bulk soils

####### Make key linking Plant_IDs to their Sites #######
smd<-as(sample_data(thresh1fieldEco_99),'data.frame') %>% select(Plant_ID,Site) %>% unique()
####### Make key linking OTUs to their taxonomic info #######
tax<-as.data.frame(as(tax_table(thresh1fieldEco_99),'matrix'))
tax$OTU_ID<-rownames(tax)

####### Look at overlap by plant- what % of phyllosphere OTUs are also detected in roots of the same plant? #######
for (i in 1:length(pairs)){ # look in one plant at a time
    if (i==1){ # initialize data frame
      pairs.df<-data.frame("Plant_ID"=character(),"sharedReadsLeaf"=numeric(),"allReadsLeaf"=numeric(),
                           "sharedReadsRoot"=numeric(),"allReadsRoot"=numeric()) # initialize data frame to be filled with plant-wise data 
      # (each plant will have data from 2 samples- root and leaf) 
    }
  plant<-pairs[i] # go through one plant at a time
  leafsubset<-subset_samples(leaf.paired,Plant_ID==plant) %>% prune_taxa(taxa_sums(.)>0,.) # look at leaf OTU table for this plant only and remove OTUs not present
  rootsubset<-subset_samples(root.paired,Plant_ID==plant) %>% prune_taxa(taxa_sums(.)>0,.) # look at root OTU table for this plant only and remove OTUs not present
  overlap<-intersect(taxa_names(leafsubset),taxa_names(rootsubset)) # identify OTUs present in both leaf and root
  leafonly<-setdiff(taxa_names(leafsubset),taxa_names(rootsubset)) # identify OTUs present in leaf only
  sharedReadsLeaf<-prune_taxa(overlap,leafsubset) %>% taxa_sums() %>% sum()
  allReadsLeaf<- taxa_sums(leafsubset) %>% sum()
  sharedReadsRoot<-prune_taxa(overlap,rootsubset) %>% taxa_sums() %>% sum()
  allReadsRoot<- taxa_sums(rootsubset) %>% sum()
  pairs.df<-rbind(pairs.df,data.frame("Plant_ID"=plant,"sharedReadsLeaf"=sharedReadsLeaf,"allReadsLeaf"=allReadsLeaf,
                                      "sharedReadsRoot"=sharedReadsRoot,"allReadsRoot"=allReadsRoot))
}

# Calculate percentages of reads in leaf samples & in root samples that are common to both organs
pairs.df<-mutate(pairs.df,
                 PctSharedLeaf=sharedReadsLeaf/allReadsLeaf,
                 PctSharedRoot=sharedReadsRoot/allReadsRoot) %>%
  merge(.,smd,by='Plant_ID')
mean(pairs.df$PctSharedLeaf) # on average, 74.5% of leaf OTUs (weighted by abundance) also found in roots
sd(pairs.df$PctSharedLeaf) # st dev = 14.5%

mean(pairs.df$PctSharedRoot) # on average, 30.3% of root OTUs (weighted by abundance) also found in leaves
sd(pairs.df$PctSharedRoot) # st dev = 14.4%

####### Figure 8a: % phyllosphere OTUs also detected in roots #######

pairs.df %>% arrange(-PctSharedLeaf) %>% 
  mutate(Plant_ID=ordered(Plant_ID,levels=Plant_ID)) %>%
  gather(key=Organ,value=PctShared,PctSharedLeaf,PctSharedRoot) %>%
  mutate(Organ=gsub("PctShared","",Organ)) -> pairs.melted
  
pdf(file='plots/Fig_8a_barplot_pct_shared.pdf',width=13,height=9)
ggplot(subset(pairs.melted,Organ=='Leaf'),aes(x=Plant_ID,y=100*PctShared))+
  geom_bar(stat='identity',fill='forest green')+
  geom_bar(data=subset(pairs.melted,Organ=='Root'),stat='identity',mapping=aes(x=Plant_ID,y=100*PctShared),fill='grey')+
  theme_classic()+
  ylab("Abundance-weighted percentage of\nOTUs shared between leaves and roots")+
  xlab("Plant")+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  theme(axis.title.x = element_text( size=30,face="bold"),
        axis.title.y  = element_text(size=30,face="bold"),
        axis.text.y=element_text(size=28,face="bold"),legend.background = element_rect(fill="gray95"),
        legend.title=element_blank(),
        legend.text=element_text(size=30,face="bold"),legend.key.height=unit(2,'lines'),legend.key.width=unit(2,'lines')) 
dev.off()

####### Look at overlap by OTU- for each OTU in the dataset, how abundant is it in leaves and roots? ####### 
otulist<-taxa_names(thresh1fieldEco_99) #5,862 OTUs total

for (i in 1:length(pairs)){ # # evaluate OTU relative abundance in leaves and roots of 1 plant at a time
  print(paste0("*** Begin analysis for plant ",i," of ",length(pairs)," ***"))
  if (i==1){ # initialize data frame to be filled with OTU-wise data 
    otu99.df<-data.frame("OTU_ID"=factor(), "Plant_ID"=factor(),
                         "leafTotal"=numeric(),"rootTotal"=numeric(),
                         "leafRelAbund"=numeric(),"rootRelAbund"=numeric(),
                         "leafCount"=numeric(),"rootCount"=numeric(),"Site"=factor(),
                         "Kingdom"=factor(),"Phylum"=factor(),"Class"=factor(),
                         "Order"=factor(),"Family"=factor())  } # end initializing data frame
  plant<-pairs[i] # store Plant_ID
  site<-subset(smd,Plant_ID==plant)$Site # store this plant's Site
  # look at just one plant at a time, store vector of counts for each OTU:
  leafCount<-subset_samples(thresh1fieldEco_99,Plant_ID==plant & Type=='leaf') %>% taxa_sums(.) %>% as.matrix()
  rootCount<-subset_samples(thresh1fieldEco_99,Plant_ID==plant & Type=='root') %>% taxa_sums(.) %>% as.matrix()
  leafTotal<-sum(leafCount) # get total Reads for this plant's leaf
  rootTotal<-sum(rootCount) # get total Reads for this plant's root
  leafRelAbund<-leafCount/leafTotal # divide by totals to get relative abundances 
  rootRelAbund<-rootCount/rootTotal # 
  this.df<-as.data.frame(cbind(leafCount,rootCount,leafRelAbund,rootRelAbund))
  this.df$OTU_ID=rownames(this.df)
  this.df<-mutate(this.df, "rootTotal"=rootTotal,'leafTotal'=leafTotal,'Plant_ID'=plant,'Site'=site) %>% 
    plyr::rename(replace=c("V1"="leafCount","V2"="rootCount","V3"="leafRelAbund","V4"="rootRelAbund")) %>%
    merge(.,tax,by='OTU_ID')
  otu99.df<-rbind(otu99.df,this.df); rm(this.df)
} # end looping through all plants

otu99.df<-mutate(otu99.df, Category=ifelse(leafCount>0&rootCount>0,"shared","unshared")) # mark each OTU as shared or unshared
otu99.df<-mutate(otu99.df,OTU_ID=factor(OTU_ID),Category=factor(Category))
save(otu99.df,file="otu99pct/OTUwise.RData")

####### Figure 8b: Scatterplot: relative abundance of shared OTUs in leaves and roots #######
shared.df<-filter(otu99.df,Category=='shared') %>%  # don't count observations from plants in which the OTU was not shared
  group_by(OTU_ID) %>%  # for each OTU...
  summarize(leafRelAbund=sum(leafCount)/sum(leafTotal),rootRelAbund=sum(rootCount)/sum(rootTotal)) %>% # calculate how common it was in leaves & roots overall
  arrange(rootRelAbund)

pdf(file="plots/Fig_8b_otuwise.pdf",width=9,height=9)
mutate(shared.df,RelAbund=rootRelAbund) %>% arrange(rootRelAbund) %>%
  within(.,OTU_ID<-as.ordered(factor(OTU_ID,levels=OTU_ID[order(rootRelAbund)]))) %>%
ggplot(.,aes(x=OTU_ID,y=log(RelAbund),color="in roots"))+
  geom_point(size=2,alpha=0.6)+
  geom_point(data=mutate(shared.df,RelAbund=leafRelAbund),aes(y=log(RelAbund),x=OTU_ID,color="in leaves"),alpha=0.6,size=2)+
  scale_y_continuous(labels=function(y) paste0(format(100*exp(y),drop0trailing=TRUE,scientific=FALSE),"%"),
                   breaks=c(log(0.0000001),log(0.000001),log(0.00001),log(0.0001),log(0.001),log(0.01),log(0.1)))+
  theme_classic()+xlab("OTU (binned at 99%)")+ylab("Relative abundance")+
  scale_color_manual(values=c("in roots"="black","in leaves"="forest green"))+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  theme(axis.title.x = element_text( size=30,face="bold"),
        axis.title.y  = element_text(size=30,face="bold"),
        axis.text.y=element_text(size=28,face="bold"),legend.background = element_rect(fill="gray95"),
        legend.title=element_blank(),
        legend.text=element_text(size=30,face="bold"),legend.key.height=unit(2,'lines'),legend.key.width=unit(2,'lines')) 
dev.off()

####### What % of shared OTUs are less common in leaves than roots? #######
dim(shared.df) # 5026 OTUs total
filter(shared.df,leafRelAbund<rootRelAbund) %>% dim() # 2526 OTUs less abundant in leaves
2526/5026 # 50.3% of OTUs more common in leaves than roots

####### mixed-effects ANCOVA #######

model<-lmer(log(leafRelAbund)~log(rootRelAbund)+(1|OTU_ID)+(1|Plant_ID),data=subset(otu99.df,Category=='shared'))
anova(model)$P # P<2.2e-16
fixef(model) # coefficient: 0.088
r2.LMM(model) # r-squared = approximately 0.6129407
r2.LMM(update(model,.~.-log(rootRelAbund))) # r-squared without root relative abundance = 0.609353
# Amount of variance in leaf relative abundance explained by root relative abundance:
0.6129407 - 0.609353
# 0.0035877

####### On average, what % of root OTUs are also found in the leaf of the same plant? #######
mean(pairs.df$PctSharedRoot) # 30.32%
sd(pairs.df$PctSharedRoot) # +/- 14.4%


####### Save image #######
save.image(paste0("otu99pct/image_",date(),".RData"))

