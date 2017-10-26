### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
#####  Linear mixed models of individual taxon counts (variance stablizing transformed) 
##### to estimate heritability in each site #
" ***** NOTE: we used a computing cluster to model each taxon in parallel. 
An executable sample of the code is provided for one taxon, but full replication of our results
will require access to a high-performance computing cluster and modification of these scripts 
to fit your own working directory paths and the syntax for your cluster's job scheduler (if not SLURM).
Alternatively, this analysis can be achieved sequentially on a single machine but will take a long time."

####### Clear workspace ########
rm(list=ls())
####### Session Info #######
sessionInfo()
"R version 3.2.3 (2015-12-10)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
[10] base     

other attached packages:
[1] VennDiagram_1.6.17         futile.logger_1.4.1        BiocParallel_1.4.3        
[4] reshape2_1.4.1             lmerTest_2.0-30            lme4_1.1-11               
[7] Matrix_1.2-4               scales_0.4.0               ggplot2_2.1.0             
[10] vegan_2.3-4                lattice_0.20-33            permute_0.9-0             
[13] doParallel_1.0.10          iterators_1.0.8            foreach_1.4.3             
[16] dplyr_0.4.3                plyr_1.8.3                 genefilter_1.52.1         
[19] Biostrings_2.38.4          XVector_0.10.0             ape_3.4                   
[22] mapdata_2.2-6              maps_3.1.0                 DESeq2_1.10.1             
[25] RcppArmadillo_0.6.600.4.0  Rcpp_0.12.3                SummarizedExperiment_1.0.2
[28] Biobase_2.30.0             GenomicRanges_1.22.4       GenomeInfoDb_1.6.3        
[31] IRanges_2.4.8              S4Vectors_0.8.11           BiocGenerics_0.16.1       
[34] phyloseq_1.14.0           

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
[49] acepack_1.3-3.3     "

####### Load source file #######
source('ecotypes_source.R')
###### Load data files ######

# load variance stablizing transformed OTU tables:
load("higher_tax_levels//vst_fieldEco_fam_CNC_unfiltered.RData") # family level
load("higher_tax_levels//vst_fieldEco_ord_CNC_unfiltered.RData") # order level
load("higher_tax_levels//vst_fieldEco_cla_CNC_unfiltered.RData") # class level
load("higher_tax_levels//vst_fieldEco_phy_CNC_unfiltered.RData") # phylum level
load("intermediate_data//phylo_leaf3_withEndog_vst.RData") # OTU level- leaves, main experiment only
leaf3.vst<-subset_samples(leaf3.withEndog.vst,Age!='endog'); rm(leaf3.withEndog.vst) # exclude Endogenous samples for modeling purposes
load("intermediate_data//phylo_root3_withEndog_vst.RData") # OTU level- roots, main experiment only
root3.vst<-subset_samples(root3.withEndog.vst,Age!='endog') ; rm(root3.withEndog.vst) # exclude Endogenous samples for modeling purposes

# load untransformed fam/.../phy tables:
## We need these to determine the rarest taxa to exclude (below)
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_fam.RData") # family level
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_ord.RData") # order level
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_cla.RData") # class level
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_phy.RData") # phylum level
load("intermediate_data/phylo_leaf3_withEndog.RData") # OTU level- leaves, main experiment only
load("intermediate_data/phylo_root3_withEndog.RData") # OTU level- roots, main experiment only
leaf3<-subset_samples(leaf3.withEndog,Age!='endog')%>% prune_taxa(taxa_sums(.)>0,.); rm(leaf3.withEndog) # exclude Endogenous samples for modeling purposes
root3<-subset_samples(root3.withEndog,Age!='endog')%>% prune_taxa(taxa_sums(.)>0,.); rm(root3.withEndog) # exclude Endogenous samples for modeling purposes

# load sample metadata
load("intermediate_data//smd_leaf3_withEndog.RData")
load("intermediate_data//smd_root3_withEndog.RData")
leaf3.smd<-subset(leaf3.smd.withEndog,Age!='endog'); rm(leaf3.smd.withEndog) # exclude Endogenous samples for modeling purposes
root3.smd<-subset(root3.smd.withEndog,Age!='endog'); rm(root3.smd.withEndog) # exclude Endogenous samples for modeling purposes

####### Save taxonomic information for each OTU #######
## Load full phyloseq object (Roots and leaves) to make sure all OTUs are represented
load("intermediate_data/Phyloseq_fieldEco_vst.RData")
fieldEco.vst %>% subset_samples(Site %in% c('Jam','Mah','Sil') & Age!='endog') -> both3.vst # leaf + root samples from main experiment in same object
taxtable<-as.data.frame(tax_table(both3.vst))

####### Create subdirectory for each taxonomic level ####### 
dir.create('heritability/otu')
dir.create('heritability/fam')
dir.create('heritability/cla')
dir.create('heritability/ord')
dir.create('heritability/phy')

####### Reduce OTU tables to include only OTUs with abundances >10% of mean ("common" OTUs) #######
# Leaves: extract "non-rare" OTUs (count totals >10% of mean)
sum(sort(taxa_sums(leaf3),decreasing=TRUE)[1])/sum(taxa_sums(leaf3)) # OTU_3 (Sphingomonas sp.) alone accounts for 36.1% of all observations
median(sort(taxa_sums(leaf3),decreasing=TRUE)) # 46
mean(taxa_sums(leaf3)) # 1693.25
mean(sort(taxa_sums(leaf3),decreasing=TRUE)[-1]) # without OTU_3, mean= 1081.96

com.leaf<-rownames(subset(data.frame("Count"=taxa_sums(leaf3)),Count>0.1*mean(taxa_sums(leaf3)))) # store "common" OTU_IDs: contains 1026 leaf OTUs
sort(taxa_sums(prune_taxa(com.leaf,leaf3))/sum(taxa_sums(leaf3)), decreasing=TRUE) # relative abundance (in entire dataset) 
sum(taxa_sums(prune_taxa(com.leaf,leaf3)))/sum(taxa_sums(leaf3)) # common set represents 98.35% of all observations

# Roots: extract "non-rare" OTUs (count totals >10% of mean)
sum(sort(taxa_sums(root3),decreasing=TRUE)[1])/sum(taxa_sums(root3)) # most abundant OTU accounts for 4.94% of all observations
median(sort(taxa_sums(root3),decreasing=TRUE)) # 780
mean(taxa_sums(root3)) # 3108.95

com.root<-rownames(subset(data.frame("Count"=taxa_sums(root3)),Count>0.1*mean(taxa_sums(root3)))) # store "common" OTU_IDs: contains 2687 root OTUs
sort(taxa_sums(prune_taxa(com.root,root3))/sum(taxa_sums(root3)), decreasing=TRUE) # relative abundance (in entire dataset)
sum(taxa_sums(prune_taxa(com.root,root3)))/sum(taxa_sums(root3)) # root common set represents 98.8% of all observations

### Prune OTU tables to 'non rare' OTUs only
leaf3.common<-prune_taxa(com.leaf,leaf3); ntaxa(leaf3.common) # 1016 remain
root3.common<-prune_taxa(com.root,root3); ntaxa(root3.common) # 2666 remain
leaf3.common.vst<-prune_taxa(com.leaf,leaf3.vst)
root3.common.vst<-prune_taxa(com.root,root3.vst)

save(leaf3.common,file="intermediate_data/phylo_leaf3_common_OTUs.RData")
save(root3.common,file="intermediate_data/phylo_root3_common_OTUs.RData")

####### Fill OTU subdirectory ######
## Make data frames with sample metadata in first columns, then Family counts in remaining columns
df.leaf.common.vst<-as.data.frame(t(as(otu_table(leaf3.common.vst),"matrix"))) %>% # extract common fam table into a data frame
  merge(leaf3.smd,.,by.y="row.names",by.x="SampleID") # combine fam table with sample metadata
df.root.common.vst<-as.data.frame(t(as(otu_table(root3.common.vst),"matrix"))) %>% # extract common fam table into a data frame
  merge(root3.smd,.,by.y="row.names",by.x="SampleID") # combine fam table with sample metadata

saveRDS(df.leaf.common.vst,file="./heritability/otu/df_leaf_vst.rds")
saveRDS(df.root.common.vst,file="./heritability/otu/df_root_vst.rds")
####### Fill Family subdirectory #######
# Pare down datasets to only the main experiment (3 sites, no endogenous): 
leaf3.fam<-subset_samples(fullEco.fam,Type=='leaf'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>% prune_taxa(taxa_sums(.)>0,.)
root3.fam<-subset_samples(fullEco.fam,Type=='root'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>% prune_taxa(taxa_sums(.)>0,.)

# Exclude rarest taxa: keep taxa whose abundances are >10% of the mean abundance
com.fam.leaf<-rownames(subset(data.frame("Count"=taxa_sums(leaf3.fam)),Count>0.1*mean(taxa_sums(leaf3.fam)))) # store "common" families: contains 69 leaf families
com.fam.root<-rownames(subset(data.frame("Count"=taxa_sums(root3.fam)),Count>0.1*mean(taxa_sums(root3.fam)))) # store "common" families: contains 110 root families

# How much data did we lose? 
sum(taxa_sums(prune_taxa(com.fam.leaf,leaf3.fam)))/sum(taxa_sums(leaf3.fam)) # common leaf set represents 99.1% of all leaf observations
sum(taxa_sums(prune_taxa(com.fam.root,root3.fam)))/sum(taxa_sums(root3.fam)) # common root set represents 99.1% of all root observations

# Prune variance-stabilizing-transformed dataset to include only the common taxa identified above
leaf3.fam.vst<-subset_samples(fieldEco.fam.vst,Type=='leaf'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>%
  prune_taxa(com.fam.leaf,.)
root3.fam.vst<-subset_samples(fieldEco.fam.vst,Type=='root'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>%
  prune_taxa(com.fam.root,.)
rm(fieldEco.fam.vst, fullEco.fam) # clean up

## Make data frames with sample metadata in first columns, then Family counts in remaining columns
df.leaf.fam.vst<-as.data.frame(t(as(otu_table(leaf3.fam.vst),"matrix"))) %>% # extract common fam table into a data frame
  merge(leaf3.smd,.,by.y="row.names",by.x="SampleID") # combine fam table with sample metadata
df.root.fam.vst<-as.data.frame(t(as(otu_table(root3.fam.vst),"matrix"))) %>% # extract common fam table into a data frame
  merge(root3.smd,.,by.y="row.names",by.x="SampleID") # combine fam table with sample metadata

saveRDS(df.leaf.fam.vst,file="./heritability/fam/df_leaf_vst.rds")
saveRDS(df.root.fam.vst,file="./heritability/fam/df_root_vst.rds")
####### Fill Order subdirectory #######
# Pare down datasets to only the main experiment (3 sites, no endogenous): 
leaf3.ord<-subset_samples(fullEco.ord,Type=='leaf'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>% prune_taxa(taxa_sums(.)>0,.)
root3.ord<-subset_samples(fullEco.ord,Type=='root'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>% prune_taxa(taxa_sums(.)>0,.)

# Exclude rarest taxa: keep taxa whose abundances are >10% of the mean abundance
com.ord.leaf<-rownames(subset(data.frame("Count"=taxa_sums(leaf3.ord)),Count>0.1*mean(taxa_sums(leaf3.ord)))) # store "common" orders: contains 38 leaf orders
com.ord.root<-rownames(subset(data.frame("Count"=taxa_sums(root3.ord)),Count>0.1*mean(taxa_sums(root3.ord)))) # store "common" orders: contains 76 root orders

# How much data did we lose? 
sum(taxa_sums(prune_taxa(com.ord.leaf,leaf3.ord)))/sum(taxa_sums(leaf3.ord)) # common leaf set represents 99.1% of all leaf observations
sum(taxa_sums(prune_taxa(com.ord.root,root3.ord)))/sum(taxa_sums(root3.ord)) # common root set represents 99.0% of all root observations

# Prune variance-stabilizing-transformed dataset to include only the common taxa identified above
leaf3.ord.vst<-subset_samples(fieldEco.ord.vst,Type=='leaf'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>%
  prune_taxa(com.ord.leaf,.)
root3.ord.vst<-subset_samples(fieldEco.ord.vst,Type=='root'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>%
  prune_taxa(com.ord.root,.)
rm(fieldEco.ord.vst, fullEco.ord) # clean up

## Make data frames with sample metadata in first columns, then ordily counts in remaining columns
df.leaf.ord.vst<-as.data.frame(t(as(otu_table(leaf3.ord.vst),"matrix"))) %>% # extract common ord table into a data frame
  merge(leaf3.smd,.,by.y="row.names",by.x="SampleID") # combine ord table with sample metadata
df.root.ord.vst<-as.data.frame(t(as(otu_table(root3.ord.vst),"matrix"))) %>% # extract common ord table into a data frame
  merge(root3.smd,.,by.y="row.names",by.x="SampleID") # combine ord table with sample metadata

saveRDS(df.leaf.ord.vst,file="./heritability/ord/df_leaf_vst.rds")
saveRDS(df.root.ord.vst,file="./heritability/ord/df_root_vst.rds")
####### Fill Class subdirectory #######
# Pare down datasets to only the main experiment (3 sites, no endogenous): 
leaf3.cla<-subset_samples(fullEco.cla,Type=='leaf'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>% prune_taxa(taxa_sums(.)>0,.)
root3.cla<-subset_samples(fullEco.cla,Type=='root'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>% prune_taxa(taxa_sums(.)>0,.)

# Exclude rarest taxa: keep taxa whose abundances are >10% of the mean abundance
com.cla.leaf<-rownames(subset(data.frame("Count"=taxa_sums(leaf3.cla)),Count>0.1*mean(taxa_sums(leaf3.cla)))) # store "common" classes: contains 27 leaf classes
com.cla.root<-rownames(subset(data.frame("Count"=taxa_sums(root3.cla)),Count>0.1*mean(taxa_sums(root3.cla)))) # store "common" classes: contains 44 root classes

# How much data did we lose? 
sum(taxa_sums(prune_taxa(com.cla.leaf,leaf3.cla)))/sum(taxa_sums(leaf3.cla)) # common leaf set represents 99.2% of all leaf observations
sum(taxa_sums(prune_taxa(com.cla.root,root3.cla)))/sum(taxa_sums(root3.cla)) # common root set represents 99.0% of all root observations

# Prune variance-stabilizing-transformed dataset to include only the common taxa identified above
leaf3.cla.vst<-subset_samples(fieldEco.cla.vst,Type=='leaf'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>%
  prune_taxa(com.cla.leaf,.)
root3.cla.vst<-subset_samples(fieldEco.cla.vst,Type=='root'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>%
  prune_taxa(com.cla.root,.)
rm(fieldEco.cla.vst, fullEco.cla) # clean up

## Make data frames with sample metadata in first columns, then claily counts in remaining columns
df.leaf.cla.vst<-as.data.frame(t(as(otu_table(leaf3.cla.vst),"matrix"))) %>% # extract common cla table into a data frame
  merge(leaf3.smd,.,by.y="row.names",by.x="SampleID") # combine cla table with sample metadata
df.root.cla.vst<-as.data.frame(t(as(otu_table(root3.cla.vst),"matrix"))) %>% # extract common cla table into a data frame
  merge(root3.smd,.,by.y="row.names",by.x="SampleID") # combine cla table with sample metadata

saveRDS(df.leaf.cla.vst,file="./heritability/cla/df_leaf_vst.rds")
saveRDS(df.root.cla.vst,file="./heritability/cla/df_root_vst.rds")
####### Fill Phylum subdirectory ########
# Pare down datasets to only the main experiment (3 sites, no endogenous): 
leaf3.phy<-subset_samples(fullEco.phy,Type=='leaf'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>% prune_taxa(taxa_sums(.)>0,.)
root3.phy<-subset_samples(fullEco.phy,Type=='root'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>% prune_taxa(taxa_sums(.)>0,.)

# Exclude rarest taxa: keep taxa whose abundances are >10% of the mean abundance
com.phy.leaf<-rownames(subset(data.frame("Count"=taxa_sums(leaf3.phy)),Count>0.1*mean(taxa_sums(leaf3.phy)))) # store "common" phyla: contains 10 leaf phyla
com.phy.root<-rownames(subset(data.frame("Count"=taxa_sums(root3.phy)),Count>0.1*mean(taxa_sums(root3.phy)))) # store "common" phyla: contains 12 root phyla

# How much data did we lose? 
sum(taxa_sums(prune_taxa(com.phy.leaf,leaf3.phy)))/sum(taxa_sums(leaf3.phy)) # common leaf set represents 99.5% of all leaf observations
sum(taxa_sums(prune_taxa(com.phy.root,root3.phy)))/sum(taxa_sums(root3.phy)) # common root set represents 98.9% of all root observations

# Prune variance-stabilizing-transformed dataset to include only the common taxa identified above
leaf3.phy.vst<-subset_samples(fieldEco.phy.vst,Type=='leaf'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>%
  prune_taxa(com.phy.leaf,.)
root3.phy.vst<-subset_samples(fieldEco.phy.vst,Type=='root'&Site%in%c('Jam','Mah','Sil')&Age!='endog') %>%
  prune_taxa(com.phy.root,.)
rm(fieldEco.phy.vst, fullEco.phy) # clean up

## Make data frames with sample metadata in first columns, then phyily counts in remaining columns
df.leaf.phy.vst<-as.data.frame(t(as(otu_table(leaf3.phy.vst),"matrix"))) %>% # extract common phy table into a data frame
  merge(leaf3.smd,.,by.y="row.names",by.x="SampleID") # combine phy table with sample metadata
df.root.phy.vst<-as.data.frame(t(as(otu_table(root3.phy.vst),"matrix"))) %>% # extract common phy table into a data frame
  merge(root3.smd,.,by.y="row.names",by.x="SampleID") # combine phy table with sample metadata

saveRDS(df.leaf.phy.vst,file="./heritability/phy/df_leaf_vst.rds")
saveRDS(df.root.phy.vst,file="./heritability/phy/df_root_vst.rds")

####### Save lists of the "major" taxa for future use #######
save(com.leaf,file="heritability/major_leaf_otu_list.RData")
save(com.fam.leaf,file="heritability/major_leaf_fam_list.RData")
save(com.ord.leaf,file="heritability/major_leaf_ord_list.RData")
save(com.cla.leaf,file="heritability/major_leaf_cla_list.RData")
save(com.phy.leaf,file="heritability/major_leaf_phy_list.RData")

save(com.root,file="heritability/major_root_otu_list.RData")
save(com.fam.root,file="heritability/major_root_fam_list.RData")
save(com.ord.root,file="heritability/major_root_ord_list.RData")
save(com.cla.root,file="heritability/major_root_cla_list.RData")
save(com.phy.root,file="heritability/major_root_phy_list.RData")

####### Estimate heritability of each taxon in each site + across all sites #######
" ***** NOTE: here I give an example of the code for a single taxon in a single organ. Due to the large number of taxa,
for replication I strongly recommend adapting this code to run in parallel on a high performance computing cluster
This script (heritability.R) includes a description of our approach, and we also provide the auxiliary scripts
that we executed on the cluster (populate_vstLMM.R, countmodels_vstLMM.R). They will need to be modified to match your own working directory paths and
the commands used by your cluster's job scheduler, if different from SLURM. Alternatively, adapt the below code
to run in sequence for all taxa on a single machine, but be advised it will take a long time to complete. "

####### This example estimates heritability of Proteobacteria abundance in leaves. #######

# Fit random-effects model for each site
lmm.rand.Jam<-lmer(Proteobacteria~(1|Genotype)+(1|Age)+(1|Harvested)+(1|Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=subset(df.leaf.phy.vst,Site=='Jam'),REML=TRUE)
lmm.rand.Mah<-lmer(Proteobacteria~(1|Genotype)+(1|Age)+(1|Harvested)+(1|Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=subset(df.leaf.phy.vst,Site=='Mah'),REML=TRUE)
lmm.rand.Sil<-lmer(Proteobacteria~(1|Genotype)+(1|Age)+(1|Harvested)+(1|Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=subset(df.leaf.phy.vst,Site=='Sil'),REML=TRUE)

lmer.rand.all<-lmer(Proteobacteria~(1|Site)+(1|Genotype)+(1|Age)+(1|Harvested)+(1|Site:Genotype)+(1|Site:Age)+(1|Site:Harvested)+(1|Site:Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=df.leaf.phy.vst,REML=TRUE)

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

# ------- * our approach to execute random-effects models in parallel on high performance computing cluster * ------- 

"
1. on local machine: 
cd ./heritability
for level in otu fam ord cla phy; do
cp ./countmodels_vstLMM.R $level/
cp ./populate_vstLMM.R $level/
done

2. cd to working directory on cluster, then:

for level in otu fam ord cla phy; do
mkdir $level
done

# transfer data files and scripts from local machine to HPCC folders for each taxonomic level:

put <this R project directory>/heritability/otu/* otu/
put <this R project directory>/heritability/fam/* fam/
put <this R project directory>/heritability/ord/* ord/
put <this R project directory>/heritability/cla/* cla/
put <this R project directory>/heritability/phy/* phy/

3. ssh on HPCC: cd <cluster working directory>
# Set up a folder for each taxon and start bash scripts for parallel job
for level in otu fam ord cla phy; do
cd $level; Rscript populate_vstLMM.R
cd ..
done

4. on cluster: 
# WHEN ALL JOBS ARE DONE: gather results:

cd <cluster working directory>
mkdir output
for level in otu fam ord cla phy; do
cd $level

# Concatenate heritability estimates from all taxa at this level (separate files for root, leaf) and store in output/
cat root/*/*vcov.txt > ../output/$level.h2.root.wHeaders.vcov.txt
head -1 ../output/$level.h2.root.wHeaders.vcov.txt > ../output/$level.h2.root.vcov.txt
grep -v Taxon* ../output/$level.h2.root.wHeaders.vcov.txt >> ../output/$level.h2.root.vcov.txt
rm ../output/$level.h2.root.wHeaders.vcov.txt

cat leaf/*/*vcov.txt > ../output/$level.h2.leaf.wHeaders.vcov.txt
head -1 ../output/$level.h2.leaf.wHeaders.vcov.txt > ../output/$level.h2.leaf.vcov.txt
grep -v Taxon* ../output/$level.h2.leaf.wHeaders.vcov.txt >> ../output/$level.h2.leaf.vcov.txt
rm ../output/$level.h2.leaf.wHeaders.vcov.txt

cd <cluster working directory>
done

5. sftp transfer concatenated results files to local machine.
cd <cluster working directory>/output
get otu*.txt <this R project directory>/heritability/otu
get fam*.txt <this R project directory>/heritability/fam
get ord*.txt <this R project directory>/heritability/ord
get cla*.txt <this R project directory>/heritability/cla
get phy*.txt <this R project directory>/heritability/phy
"

# Back to R:

####### Load heritability estimates and plot for each site & across sites #######
h2.vcov.taxa<-read.delim('heritability/otu/otu.h2.leaf.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='leaf',Level='OTU')
read.delim('heritability/fam/fam.h2.leaf.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='leaf',Level='Family') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
read.delim('heritability/ord/ord.h2.leaf.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='leaf',Level='Order') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
read.delim('heritability/cla/cla.h2.leaf.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='leaf',Level='Class') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
read.delim('heritability/phy/phy.h2.leaf.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='leaf',Level='Phylum') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
read.delim('heritability/otu/otu.h2.root.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='root',Level='OTU') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
read.delim('heritability/fam/fam.h2.root.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='root',Level='Family') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
read.delim('heritability/ord/ord.h2.root.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='root',Level='Order') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
read.delim('heritability/cla/cla.h2.root.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='root',Level='Class') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
read.delim('heritability/phy/phy.h2.root.vcov.txt',sep='\t',header=TRUE) %>% mutate(Tissue='root',Level='Phylum') %>% rbind(.,h2.vcov.taxa) -> h2.vcov.taxa
h2.vcov.taxa<-rename(h2.vcov.taxa,H2=h2_vcov)
# combine with phylogenetic information
h2.vcov.otus<-merge(subset(h2.vcov.taxa,Level=='OTU'),taxtable,by.x="Taxon",by.y="row.names",suffixes=c("",".y"))

####### Fig. 6a: Boxplot: heritability estimates #######
# Just show the major phyla
mainphyla<-c("Proteobacteria","Actinobacteria","Acidobacteria","Bacteroidetes","Cyanobacteria","Crenarchaeota","Chloroflexi","Armatimonadetes","Verrucomicrobia","Planctomycetes","Gemmatimonadetes","Firmicutes")

pdf(file="plots/Fig_6a_boxplot_heritability.pdf",width=15,height=9)
update_geom_defaults("point", list(colour = NULL))
ggplot(data = subset(h2.vcov.otus,Phylum%in%mainphyla),aes(x=Phylum,y=H2,colour=Site))+
  geom_boxplot(size=0.8,outlier.colour=NULL,position=position_dodge(w=0.7),outlier.size=2)+
  facet_grid(Tissue~.)+
  ylab("Broad-sense heritability (%)\n")+
  scale_color_manual(values=c("red",sitePalette),name="Site:   ")+
  theme_minimal()+scale_y_continuous(labels=function(y) 100*y)+
  theme(axis.text.x = element_text(angle=25,vjust=1, hjust=0.95,size=28,face="bold"))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=34,face="bold"))+
  theme(axis.text.y=element_text(size=28,face="bold"))+
  theme(panel.margin=unit(2,"lines"),legend.position="top")+
  theme(strip.background=element_rect(color="gray90",fill="gray90"))+
  theme(legend.title=element_text(size=36,face="bold"),legend.text=element_text(size=32,face="bold"))+
  theme(legend.key.height=unit(2,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(strip.text.y = element_text(size=36,face="bold"))
update_geom_defaults("point", list(colour = "black"))
dev.off()

####### Save image #######
save.image(paste0("heritability/image_",date(),".RData"))
