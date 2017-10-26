### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
#####  Negative binomial models to estimate fold changes of individual taxa
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
####### Load source file #######
source('ecotypes_source.R')
###### Load data files ######
# load untransformed Phyloseq objects:
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_fam.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_ord.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_cla.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_phy.RData")
load("intermediate_data/phylo_leaf3_withEndog.RData")
load("intermediate_data/phylo_root3_withEndog.RData")

####### Pare down phyloseq objects to include only relevant samples + 'major' taxa (see Heritability.R) #######
# Load lists of major taxa at all levels:
load("heritability/major_leaf_otu_list.RData")
load("heritability/major_leaf_fam_list.RData")
load("heritability/major_leaf_ord_list.RData")
load("heritability/major_leaf_cla_list.RData")
load("heritability/major_leaf_phy_list.RData")

load("heritability/major_root_otu_list.RData")
load("heritability/major_root_fam_list.RData")
load("heritability/major_root_ord_list.RData")
load("heritability/major_root_cla_list.RData")
load("heritability/major_root_phy_list.RData")

# Prune out extraneous samples and minor taxa: 
leaf3.common<-subset_samples(leaf3.withEndog,Age!='endog') %>% prune_taxa(com.leaf,.)
leaf3.common.fam<-subset_samples(fullEco.fam,Site%in%c('Jam','Mah','Sil')&Type=='leaf'&Age!='endog') %>% prune_taxa(com.fam.leaf,.)
leaf3.common.ord<-subset_samples(fullEco.ord,Site%in%c('Jam','Mah','Sil')&Type=='leaf'&Age!='endog') %>% prune_taxa(com.ord.leaf,.)
leaf3.common.cla<-subset_samples(fullEco.cla,Site%in%c('Jam','Mah','Sil')&Type=='leaf'&Age!='endog') %>% prune_taxa(com.cla.leaf,.)
leaf3.common.phy<-subset_samples(fullEco.phy,Site%in%c('Jam','Mah','Sil')&Type=='leaf'&Age!='endog') %>% prune_taxa(com.phy.leaf,.)

root3.common<-subset_samples(root3.withEndog,Age!='endog') %>% prune_taxa(com.root,.)
root3.common.fam<-subset_samples(fullEco.fam,Site%in%c('Jam','Mah','Sil')&Type=='root'&Age!='endog') %>% prune_taxa(com.fam.root,.)
root3.common.ord<-subset_samples(fullEco.ord,Site%in%c('Jam','Mah','Sil')&Type=='root'&Age!='endog') %>% prune_taxa(com.ord.root,.)
root3.common.cla<-subset_samples(fullEco.cla,Site%in%c('Jam','Mah','Sil')&Type=='root'&Age!='endog') %>% prune_taxa(com.cla.root,.)
root3.common.phy<-subset_samples(fullEco.phy,Site%in%c('Jam','Mah','Sil')&Type=='root'&Age!='endog') %>% prune_taxa(com.phy.root,.)

rm(leaf3.withEndog,root3.withEndog,fullEco.fam,fullEco.ord,fullEco.cla,fullEco.phy)

sample_data(leaf3.common)<-sample_data(leaf3.common.fam) # make sample_data the same at all levels
sample_data(root3.common)<-sample_data(root3.common.fam) # make sample_data the same at all levels

#### ------ *** Begin DESeq2 log 2 fold change estimates *** ------ ####

## Negative binomial GLMs of "common" OTUs only (>10% of mean abundance)
## Justification: throwing out "rare" OTUs in that step reduced total # of observations 
## by only 1%. 

system("mkdir foldchange/DESeq2_datasets") # make directory to store DESeq2 analyses
register(MulticoreParam(),default=TRUE)
multicoreWorkers() # how many cores available?

####### All levels: Create DESeq2 datasets for MAIN EFFECTS (no interaction terms) #######
## See DESeq2 vignette for explanation of why including interaction terms alters interpretation of main effects (justification for doing 2 separate models: main effects only and interactions only)

## OTU level: roots then leaves: 
dds.root.main.otu<-phyloseq_to_deseq2(root3.common,~Site+Genotype+Age+Harvested)
dds.root.main.otu = estimateSizeFactors(dds.root.main.otu, geoMeans = apply(counts(dds.root.main.otu), 1, gm_mean))
dds.root.main.otu<-DESeq(dds.root.main.otu,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

dds.leaf.main.otu<-phyloseq_to_deseq2(leaf3.common,~Site+Genotype+Age+Harvested)
dds.leaf.main.otu = estimateSizeFactors(dds.leaf.main.otu, geoMeans = apply(counts(dds.leaf.main.otu), 1, gm_mean))
dds.leaf.main.otu<-DESeq(dds.leaf.main.otu,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

## Family level: roots then leaves: 
dds.root.main.fam<-phyloseq_to_deseq2(root3.common.fam,~Site+Genotype+Age+Harvested)
dds.root.main.fam = estimateSizeFactors(dds.root.main.fam, geoMeans = apply(counts(dds.root.main.fam), 1, gm_mean))
dds.root.main.fam<-DESeq(dds.root.main.fam,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

dds.leaf.main.fam<-phyloseq_to_deseq2(leaf3.common.fam,~Site+Genotype+Age+Harvested)
dds.leaf.main.fam = estimateSizeFactors(dds.leaf.main.fam, geoMeans = apply(counts(dds.leaf.main.fam), 1, gm_mean))
dds.leaf.main.fam<-DESeq(dds.leaf.main.fam,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

## Order level: roots then leaves: 
dds.root.main.ord<-phyloseq_to_deseq2(root3.common.ord,~Site+Genotype+Age+Harvested)
dds.root.main.ord = estimateSizeFactors(dds.root.main.ord, geoMeans = apply(counts(dds.root.main.ord), 1, gm_mean))
dds.root.main.ord<-DESeq(dds.root.main.ord,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

dds.leaf.main.ord<-phyloseq_to_deseq2(leaf3.common.ord,~Site+Genotype+Age+Harvested)
dds.leaf.main.ord = estimateSizeFactors(dds.leaf.main.ord, geoMeans = apply(counts(dds.leaf.main.ord), 1, gm_mean))
dds.leaf.main.ord<-DESeq(dds.leaf.main.ord,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

## Class level: roots then leaves: 
dds.root.main.cla<-phyloseq_to_deseq2(root3.common.cla,~Site+Genotype+Age+Harvested)
dds.root.main.cla = estimateSizeFactors(dds.root.main.cla, geoMeans = apply(counts(dds.root.main.cla), 1, gm_mean))
dds.root.main.cla<-DESeq(dds.root.main.cla,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

dds.leaf.main.cla<-phyloseq_to_deseq2(leaf3.common.cla,~Site+Genotype+Age+Harvested)
dds.leaf.main.cla = estimateSizeFactors(dds.leaf.main.cla, geoMeans = apply(counts(dds.leaf.main.cla), 1, gm_mean))
dds.leaf.main.cla<-DESeq(dds.leaf.main.cla,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

## Phylum level: roots then leaves: 
dds.root.main.phy<-phyloseq_to_deseq2(root3.common.phy,~Site+Genotype+Age+Harvested)
dds.root.main.phy = estimateSizeFactors(dds.root.main.phy, geoMeans = apply(counts(dds.root.main.phy), 1, gm_mean))
dds.root.main.phy<-DESeq(dds.root.main.phy,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))

dds.leaf.main.phy<-phyloseq_to_deseq2(leaf3.common.phy,~Site+Genotype+Age+Harvested)
dds.leaf.main.phy = estimateSizeFactors(dds.leaf.main.phy, geoMeans = apply(counts(dds.leaf.main.phy), 1, gm_mean))
dds.leaf.main.phy<-DESeq(dds.leaf.main.phy,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5))
####### All levels: save main-effects DESeq objects #######
save(dds.leaf.main.otu,file="foldchange/DESeq2_datasets/dds_leaf_otu_main.RData")
save(dds.leaf.main.fam,file="foldchange/DESeq2_datasets/dds_leaf_fam_main.RData")
save(dds.leaf.main.ord,file="foldchange/DESeq2_datasets/dds_leaf_ord_main.RData")
save(dds.leaf.main.cla,file="foldchange/DESeq2_datasets/dds_leaf_cla_main.RData")
save(dds.leaf.main.phy,file="foldchange/DESeq2_datasets/dds_leaf_phy_main.RData")

save(dds.root.main.otu,file="foldchange/DESeq2_datasets/dds_root_otu_main.RData")
save(dds.root.main.fam,file="foldchange/DESeq2_datasets/dds_root_fam_main.RData")
save(dds.root.main.ord,file="foldchange/DESeq2_datasets/dds_root_ord_main.RData")
save(dds.root.main.cla,file="foldchange/DESeq2_datasets/dds_root_cla_main.RData")
save(dds.root.main.phy,file="foldchange/DESeq2_datasets/dds_root_phy_main.RData")

####### OTU level: likelihood ratio tests for main effects #######
# LRT for Genotype: 
nbinomLRT(dds.root.main.otu,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root') -> lrt.otu
nbinomLRT(dds.leaf.main.otu,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf')  %>% rbind(.,lrt.otu) -> lrt.otu

# LRT for Site: 
nbinomLRT(dds.root.main.otu,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root')  %>% rbind(.,lrt.otu) -> lrt.otu
nbinomLRT(dds.leaf.main.otu,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf')  %>% rbind(.,lrt.otu) -> lrt.otu

# LRT for Age: 
nbinomLRT(dds.root.main.otu,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root')  %>% rbind(.,lrt.otu) -> lrt.otu
nbinomLRT(dds.leaf.main.otu,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf')  %>% rbind(.,lrt.otu) -> lrt.otu

# LRT for Year: 
nbinomLRT(dds.root.main.otu,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='root')  %>% rbind(.,lrt.otu) -> lrt.otu
nbinomLRT(dds.leaf.main.otu,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='leaf')  %>% rbind(.,lrt.otu) -> lrt.otu
####### OTU level: likelihood ratio tests for interaction effects #######
## LRT for Genotype:Site
nbinomLRT(dds.root.main.otu,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.root.main.otu)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='root')  %>% rbind(.,lrt.otu) -> lrt.otu
nbinomLRT(dds.leaf.main.otu,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.leaf.main.otu)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='leaf')  %>% rbind(.,lrt.otu) -> lrt.otu

## LRT for Age:Site
nbinomLRT(dds.root.main.otu,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.root.main.otu)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='root')  %>% rbind(.,lrt.otu) -> lrt.otu
nbinomLRT(dds.leaf.main.otu,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.leaf.main.otu)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='leaf')  %>% rbind(.,lrt.otu) -> lrt.otu

## LRT for Year:Site
nbinomLRT(dds.root.main.otu,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.root.main.otu)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='root')  %>% rbind(.,lrt.otu) -> lrt.otu
nbinomLRT(dds.leaf.main.otu,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.leaf.main.otu)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='leaf')  %>% rbind(.,lrt.otu) -> lrt.otu



















































# Analyses for entire experiment: 

####### Family level: likelihood ratio tests for main effects #######
# LRT for Genotype: 
nbinomLRT(dds.root.main.fam,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root') -> lrt.fam
nbinomLRT(dds.leaf.main.fam,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf')  %>% rbind(.,lrt.fam) -> lrt.fam

# LRT for Site: 
nbinomLRT(dds.root.main.fam,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root')  %>% rbind(.,lrt.fam) -> lrt.fam
nbinomLRT(dds.leaf.main.fam,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf')  %>% rbind(.,lrt.fam) -> lrt.fam

# LRT for Age: 
nbinomLRT(dds.root.main.fam,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root')  %>% rbind(.,lrt.fam) -> lrt.fam
nbinomLRT(dds.leaf.main.fam,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf')  %>% rbind(.,lrt.fam) -> lrt.fam

# LRT for Year: 
nbinomLRT(dds.root.main.fam,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='root')  %>% rbind(.,lrt.fam) -> lrt.fam
nbinomLRT(dds.leaf.main.fam,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='leaf')  %>% rbind(.,lrt.fam) -> lrt.fam
####### Family level: likelihood ratio tests for interaction effects #######
## LRT for Genotype:Site
nbinomLRT(dds.root.main.fam,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.root.main.fam)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='root')  %>% rbind(.,lrt.fam) -> lrt.fam
nbinomLRT(dds.leaf.main.fam,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.leaf.main.fam)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='leaf')  %>% rbind(.,lrt.fam) -> lrt.fam

## LRT for Age:Site
nbinomLRT(dds.root.main.fam,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.root.main.fam)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='root')  %>% rbind(.,lrt.fam) -> lrt.fam
nbinomLRT(dds.leaf.main.fam,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.leaf.main.fam)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='leaf')  %>% rbind(.,lrt.fam) -> lrt.fam

## LRT for Year:Site
nbinomLRT(dds.root.main.fam,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.root.main.fam)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='root')  %>% rbind(.,lrt.fam) -> lrt.fam
nbinomLRT(dds.leaf.main.fam,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.leaf.main.fam)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='leaf')  %>% rbind(.,lrt.fam) -> lrt.fam

####### Order level: likelihood ratio tests for main effects #######
# LRT for Genotype: 
nbinomLRT(dds.root.main.ord,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root') -> lrt.ord
nbinomLRT(dds.leaf.main.ord,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf')  %>% rbind(.,lrt.ord) -> lrt.ord

# LRT for Site: 
nbinomLRT(dds.root.main.ord,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root')  %>% rbind(.,lrt.ord) -> lrt.ord
nbinomLRT(dds.leaf.main.ord,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf')  %>% rbind(.,lrt.ord) -> lrt.ord

# LRT for Age: 
nbinomLRT(dds.root.main.ord,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root')  %>% rbind(.,lrt.ord) -> lrt.ord
nbinomLRT(dds.leaf.main.ord,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf')  %>% rbind(.,lrt.ord) -> lrt.ord

# LRT for Year: 
nbinomLRT(dds.root.main.ord,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='root')  %>% rbind(.,lrt.ord) -> lrt.ord
nbinomLRT(dds.leaf.main.ord,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='leaf')  %>% rbind(.,lrt.ord) -> lrt.ord
####### Order level: likelihood ratio tests for interaction effects #######
## LRT for Genotype:Site
nbinomLRT(dds.root.main.ord,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.root.main.ord)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='root')  %>% rbind(.,lrt.ord) -> lrt.ord
nbinomLRT(dds.leaf.main.ord,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.leaf.main.ord)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='leaf')  %>% rbind(.,lrt.ord) -> lrt.ord

## LRT for Age:Site
nbinomLRT(dds.root.main.ord,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.root.main.ord)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='root')  %>% rbind(.,lrt.ord) -> lrt.ord
nbinomLRT(dds.leaf.main.ord,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.leaf.main.ord)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='leaf')  %>% rbind(.,lrt.ord) -> lrt.ord

## LRT for Year:Site
nbinomLRT(dds.root.main.ord,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.root.main.ord)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='root')  %>% rbind(.,lrt.ord) -> lrt.ord
nbinomLRT(dds.leaf.main.ord,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.leaf.main.ord)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='leaf')  %>% rbind(.,lrt.ord) -> lrt.ord

####### Class level: likelihood ratio tests for main effects #######
# LRT for Genotype: 
nbinomLRT(dds.root.main.cla,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root') -> lrt.cla
nbinomLRT(dds.leaf.main.cla,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf')  %>% rbind(.,lrt.cla) -> lrt.cla

# LRT for Site: 
nbinomLRT(dds.root.main.cla,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root')  %>% rbind(.,lrt.cla) -> lrt.cla
nbinomLRT(dds.leaf.main.cla,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf')  %>% rbind(.,lrt.cla) -> lrt.cla

# LRT for Age: 
nbinomLRT(dds.root.main.cla,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root')  %>% rbind(.,lrt.cla) -> lrt.cla
nbinomLRT(dds.leaf.main.cla,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf')  %>% rbind(.,lrt.cla) -> lrt.cla

# LRT for Year: 
nbinomLRT(dds.root.main.cla,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='root')  %>% rbind(.,lrt.cla) -> lrt.cla
nbinomLRT(dds.leaf.main.cla,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='leaf')  %>% rbind(.,lrt.cla) -> lrt.cla
####### Class level: likelihood ratio tests for interaction effects #######
## LRT for Genotype:Site
nbinomLRT(dds.root.main.cla,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.root.main.cla)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='root')  %>% rbind(.,lrt.cla) -> lrt.cla
nbinomLRT(dds.leaf.main.cla,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.leaf.main.cla)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='leaf')  %>% rbind(.,lrt.cla) -> lrt.cla

## LRT for Age:Site
nbinomLRT(dds.root.main.cla,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.root.main.cla)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='root')  %>% rbind(.,lrt.cla) -> lrt.cla
nbinomLRT(dds.leaf.main.cla,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.leaf.main.cla)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='leaf')  %>% rbind(.,lrt.cla) -> lrt.cla

## LRT for Year:Site
nbinomLRT(dds.root.main.cla,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.root.main.cla)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='root')  %>% rbind(.,lrt.cla) -> lrt.cla
nbinomLRT(dds.leaf.main.cla,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.leaf.main.cla)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='leaf')  %>% rbind(.,lrt.cla) -> lrt.cla

####### Phylum level: likelihood ratio tests for main effects #######
# LRT for Genotype: 
nbinomLRT(dds.root.main.phy,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root') -> lrt.phy
nbinomLRT(dds.leaf.main.phy,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf')  %>% rbind(.,lrt.phy) -> lrt.phy

# LRT for Site: 
nbinomLRT(dds.root.main.phy,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root')  %>% rbind(.,lrt.phy) -> lrt.phy
nbinomLRT(dds.leaf.main.phy,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf')  %>% rbind(.,lrt.phy) -> lrt.phy

# LRT for Age: 
nbinomLRT(dds.root.main.phy,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root')  %>% rbind(.,lrt.phy) -> lrt.phy
nbinomLRT(dds.leaf.main.phy,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf')  %>% rbind(.,lrt.phy) -> lrt.phy

# LRT for Year: 
nbinomLRT(dds.root.main.phy,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='root')  %>% rbind(.,lrt.phy) -> lrt.phy
nbinomLRT(dds.leaf.main.phy,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='leaf')  %>% rbind(.,lrt.phy) -> lrt.phy
####### Phylum level: likelihood ratio tests for interaction effects #######
## LRT for Genotype:Site
nbinomLRT(dds.root.main.phy,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.root.main.phy)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='root')  %>% rbind(.,lrt.phy) -> lrt.phy
nbinomLRT(dds.leaf.main.phy,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.leaf.main.phy)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='leaf')  %>% rbind(.,lrt.phy) -> lrt.phy

## LRT for Age:Site
nbinomLRT(dds.root.main.phy,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.root.main.phy)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='root')  %>% rbind(.,lrt.phy) -> lrt.phy
nbinomLRT(dds.leaf.main.phy,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.leaf.main.phy)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='leaf')  %>% rbind(.,lrt.phy) -> lrt.phy

## LRT for Year:Site
nbinomLRT(dds.root.main.phy,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.root.main.phy)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='root')  %>% rbind(.,lrt.phy) -> lrt.phy
nbinomLRT(dds.leaf.main.phy,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.leaf.main.phy)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='leaf')  %>% rbind(.,lrt.phy) -> lrt.phy

####### Combine & save likelihood ratio test results for all levels, terms #######
lrt.all<-rbind(mutate(lrt.otu,Level='otu'),
               mutate(lrt.fam,Level='fam'),
               mutate(lrt.ord,Level='ord'),
               mutate(lrt.cla,Level='cla'),
               mutate(lrt.phy,Level='phy')) %>%
  rename(Taxon=OTU_ID) %>%
  mutate(Taxon=factor(Taxon),Term=factor(Term),Organ=factor(Organ),Level=factor(Level))
save(lrt.all,file="foldchange/LRT_results_all.RData")

####### All levels: Calculate log2 fold changes due to main effects #######
l2fc<-data.frame("baseMean"=numeric(),"log2FoldChange"=numeric(),"lfcSE"=numeric(),"stat"=numeric(),"pvalue"=numeric(),
           "padj"=numeric(),"Taxon"=factor(),"Term"=factor(),"Organ"=factor(),"Level"=factor())  # initalize log2 fold change results

## Get all pairwise Genotype contrasts for roots then leaves
GenoList<-c('JAM','MAH','MIL','PAR','SIL')
LevelList<-c('otu','fam','ord','cla','phy')
for (level in LevelList) {
  leafdds<-get(paste0('dds.leaf.main.',level))
  rootdds<-get(paste0('dds.root.main.',level))
  for (j in 1:(length(GenoList)-1)) { # j = index of 1st genotype
    geno1<-GenoList[j]
      for (k in (j+1):length(GenoList)) { # k = index of 2nd genotype
        geno2<-GenoList[k]
        # print(paste0(geno1,".",geno2)) # debugging line to make sure all pairwise comparisons are represented
        # get log2 fold change between this pair of genotypes
        ## roots first:
        results(rootdds,contrast=c('Genotype',geno1,geno2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
          mutate(Taxon=row.names(.),Term='Genotype',Contrast=paste0(geno1,'_',geno2),Site='all',Organ='root',Level=level) %>%
          rbind(.,l2fc) -> l2fc
        ## then leaves:
        results(leafdds,contrast=c('Genotype',geno1,geno2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
          mutate(Taxon=row.names(.),Term='Genotype',Contrast=paste0(geno1,'_',geno2),Site='all',Organ='leaf',Level=level) %>%
          rbind(.,l2fc) -> l2fc
      }
  }
  rm(leafdds,rootdds) # clean up
} 
rm(geno1,geno2,level,j,k) # clean up

## Get all pairwise Site contrasts for roots then leaves
SiteList<-c('Jam','Mah','Sil')
LevelList<-c('otu','fam','ord','cla','phy')
for (level in LevelList) {
  leafdds<-get(paste0('dds.leaf.main.',level))
  rootdds<-get(paste0('dds.root.main.',level))
  for (j in 1:(length(SiteList)-1)) { # j = index of 1st site
    site1<-SiteList[j]
    for (k in (j+1):length(SiteList)) { # k = index of 2nd site
      site2<-SiteList[k]
      # print(paste0(site1,".",site2)) # debugging line to make sure all pairwise comparisons are represented
      # get log2 fold change between this pair of sites
      ## roots first:
      results(rootdds,contrast=c('Site',site1,site2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
        mutate(Taxon=row.names(.),Term='Site',Contrast=paste0(site1,'_',site2),Site='all',Organ='root',Level=level) %>%
        rbind(.,l2fc) -> l2fc
      ## then leaves:
      results(leafdds,contrast=c('Site',site1,site2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
        mutate(Taxon=row.names(.),Term='Site',Contrast=paste0(site1,'_',site2),Site='all',Organ='leaf',Level=level) %>%
        rbind(.,l2fc) -> l2fc
    }
  }
  rm(leafdds,rootdds) # clean up
} 
rm(site1,site2,j,k,level) # clean up


# Get all pairwise Age contrasts:
AgeList<-c('2','3','4')
LevelList<-c('otu','fam','ord','cla','phy')
for (level in LevelList) {
  leafdds<-get(paste0('dds.leaf.main.',level))
  rootdds<-get(paste0('dds.root.main.',level))
  for (j in 1:(length(AgeList)-1)) { # j = index of 1st Age
    Age1<-AgeList[j]
    for (k in (j+1):length(AgeList)) { # k = index of 2nd Age
      Age2<-AgeList[k]
      # print(paste0(Age1,".",Age2)) # debugging line to make sure all pairwise comparisons are represented
      # get log2 fold change between this pair of Ages
      ## roots first:
      results(rootdds,contrast=c('Age',Age1,Age2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
        mutate(Taxon=row.names(.),Term='Age',Contrast=paste0(Age1,'_',Age2),Site='all',Organ='root',Level=level) %>%
        rbind(.,l2fc) -> l2fc
      ## then leaves:
      results(leafdds,contrast=c('Age',Age1,Age2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
        mutate(Taxon=row.names(.),Term='Age',Contrast=paste0(Age1,'_',Age2),Site='all',Organ='leaf',Level=level) %>%
        rbind(.,l2fc) -> l2fc
    }
  }
  rm(leafdds,rootdds) # clean up
} 
rm(Age1,Age2,j,k,level) # clean up


# Get Year effect for roots and leaves:
LevelList<-c('otu','fam','ord','cla','phy')
for (level in LevelList) {
  leafdds<-get(paste0('dds.leaf.main.',level))
  rootdds<-get(paste0('dds.root.main.',level))
  # get log2 fold change between Years
  ## roots first:
  results(rootdds,contrast=c('Harvested','2011','2012'),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
    mutate(Taxon=row.names(.),Term='Year',Contrast='2011_2012',Site='all',Organ='root',Level=level) %>%
    rbind(.,l2fc) -> l2fc
  ## then leaves:
  results(leafdds,contrast=c('Harvested','2011','2012'),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
    mutate(Taxon=row.names(.),Term='Year',Contrast='2011_2012',Site='all',Organ='leaf',Level=level) %>%
    rbind(.,l2fc) -> l2fc
  rm(leafdds,rootdds) # clean up
    }
rm(level) # clean up

####### All levels: site-specific log2 fold changes (interaction terms) #######
l2fc.int<-data.frame("baseMean"=numeric(),"log2FoldChange"=numeric(),"lfcSE"=numeric(),"stat"=numeric(),"pvalue"=numeric(),
                 "padj"=numeric(),"Taxon"=factor(),"Term"=factor(),"Organ"=factor(),"Level"=factor())  # initalize log2 fold change results
# rename common OTU phyloseq objects for convenience:
leaf3.common.otu<-leaf3.common ; rm(leaf3.common)
root3.common.otu<-root3.common; rm(root3.common)
####### All levels: log2 fold changes: Genotype x Site #######

LevelList<-c('otu','fam','ord','cla','phy')
GenoList<-c('JAM','MAH','MIL','PAR','SIL')
SiteList<-c('Jam','Mah','Sil')
for (level in LevelList) {
  ## Make new DESeq2 objects to look at GxS while controlling for Age, Year: 
  print(level)
  leafphylo<-get(paste0('leaf3.common.',level))
  rootphylo<-get(paste0('root3.common.',level))
  leafGxS<-phyloseq_to_deseq2(leafphylo,~1) # make DESeq2 object from phyloseq object (common taxa)
  rootGxS<-phyloseq_to_deseq2(rootphylo,~1) # make DESeq2 object from phyloseq object (common taxa)
  leafGxS$GxS <- factor(paste0(leafGxS$Genotype,".",leafGxS$Site)) # Make separate "grouping factor" with Genotype and Site
  rootGxS$GxS <- factor(paste0(rootGxS$Genotype,".",rootGxS$Site)) # Make separate "grouping factor" with Genotype and Site
  design(leafGxS) <- ~GxS+Age+Harvested # define experimental design: interested in Geno x Site averaged over ages, years
  design(rootGxS) <- ~GxS+Age+Harvested # define experimental design: interested in Geno x Site averaged over ages, years
  leafGxS = estimateSizeFactors(leafGxS, geoMeans = apply(counts(leafGxS), 1, gm_mean)) # Estimate size factors
  rootGxS = estimateSizeFactors(rootGxS, geoMeans = apply(counts(rootGxS), 1, gm_mean)) # Estimate size factors
  leafGxS<-DESeq(leafGxS,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) # Store differential abundance analysis for this level
  rootGxS<-DESeq(rootGxS,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) # Store differential abundance analysis for this level
  save(leafGxS,file=paste0("foldchange/DESeq2_datasets/dds_leaf_",level,"_GxS.RData")) # save to file
  save(rootGxS,file=paste0("foldchange/DESeq2_datasets/dds_root_",level,"_GxS.RData")) # save to file
  ## Get all pairwise Genotype contrasts  at each site
  for (site in SiteList) { 
    for (j in 1:(length(GenoList)-1)) { # j = index of 1st genotype
      geno1<-paste0(GenoList[j],".",site)
      for (k in (j+1):length(GenoList)) { # k = index of 2nd genotype
        geno2<-paste0(GenoList[k],".",site)
        # print(paste0(geno1,".",geno2)) # debugging line to make sure all pairwise comparisons are represented at each site
        # get site-specific log2 fold change between this pair of genotypes at this site
        ## roots first:
        results(rootGxS,contrast=c('GxS',geno1,geno2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
          mutate(Taxon=row.names(.),Term='GxS',Contrast=paste0(geno1,"_",geno2),Site=site,Organ='root',Level=level) %>%
          rbind(.,l2fc.int) -> l2fc.int
        ## then leaves:
        results(leafGxS,contrast=c('GxS',geno1,geno2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
          mutate(Taxon=row.names(.),Term='GxS',Contrast=paste0(geno1,"_",geno2),Site=site,Organ='leaf',Level=level) %>%
          rbind(.,l2fc.int) -> l2fc.int
      }
    }
  }
}
rm(site,level,geno1,geno2,j,k,leafGxS,rootGxS,leafphylo,rootphylo)
####### All levels: log2 fold changes: Age x Site #######

LevelList<-c('otu','fam','ord','cla','phy')
AgeList<-c('2','3','4')
SiteList<-c('Jam','Mah','Sil')
for (level in LevelList) {
  ## Make new DESeq2 objects to look at AxS while controlling for Genotype, Year: 
  print(level)
  leafphylo<-get(paste0('leaf3.common.',level))
  rootphylo<-get(paste0('root3.common.',level))
  leafAxS<-phyloseq_to_deseq2(leafphylo,~1) # make DESeq2 object from phyloseq object (common taxa)
  rootAxS<-phyloseq_to_deseq2(rootphylo,~1) # make DESeq2 object from phyloseq object (common taxa)
  leafAxS$AxS <- factor(paste0(leafAxS$Age,".",leafAxS$Site)) # Make separate "grouping factor" with Age and Site
  rootAxS$AxS <- factor(paste0(rootAxS$Age,".",rootAxS$Site)) # Make separate "grouping factor" with Age and Site
  design(leafAxS) <- ~AxS+Genotype+Harvested # define experimental design: interested in Age x Site averaged over genotypes, years
  design(rootAxS) <- ~AxS+Genotype+Harvested # define experimental design: interested in Age x Site averaged over genotypes, years
  leafAxS = estimateSizeFactors(leafAxS, geoMeans = apply(counts(leafAxS), 1, gm_mean)) # Estimate size factors
  rootAxS = estimateSizeFactors(rootAxS, geoMeans = apply(counts(rootAxS), 1, gm_mean)) # Estimate size factors
  leafAxS<-DESeq(leafAxS,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) # Store differential abundance analysis for this level
  rootAxS<-DESeq(rootAxS,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) # Store differential abundance analysis for this level
  save(leafAxS,file=paste0("foldchange/DESeq2_datasets/dds_leaf_",level,"_AxS.RData")) # save to file
  save(rootAxS,file=paste0("foldchange/DESeq2_datasets/dds_root_",level,"_AxS.RData")) # save to file
  ## Get all pairwise Age contrasts  at each site
  for (site in SiteList) { 
    for (j in 1:(length(AgeList)-1)) { # j = index of 1st age
      grp1<-paste0(AgeList[j],".",site)
      for (k in (j+1):length(AgeList)) { # k = index of 2nd age
        grp2<-paste0(AgeList[k],".",site)
        # print(paste0(grp1,".",grp2)) # debugging line to make sure all pairwise comparisons are represented at each site
        # get site-specific log2 fold change between this pair of age groups  at this site
        ## roots first:
        results(rootAxS,contrast=c('AxS',grp1,grp2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
          mutate(Taxon=row.names(.),Term='AxS',Contrast=paste0(grp1,"_",grp2),Site=site,Organ='root',Level=level) %>%
          rbind(.,l2fc.int) -> l2fc.int
        ## then leaves:
        results(leafAxS,contrast=c('AxS',grp1,grp2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
          mutate(Taxon=row.names(.),Term='AxS',Contrast=paste0(grp1,"_",grp2),Site=site,Organ='leaf',Level=level) %>%
          rbind(.,l2fc.int) -> l2fc.int
      }
    }
  }
}
rm(site,level,grp1,grp2,j,k,leafAxS,rootAxS,leafphylo,rootphylo)
####### All levels: log2 fold changes: Year x Site #######

LevelList<-c('otu','fam','ord','cla','phy')
SiteList<-c('Jam','Mah','Sil')
for (level in LevelList) {
  ## Make new DESeq2 objects to look at YxS while controlling for Genotype, Age: 
  print(level)
  leafphylo<-get(paste0('leaf3.common.',level))
  rootphylo<-get(paste0('root3.common.',level))
  leafYxS<-phyloseq_to_deseq2(leafphylo,~1) # make DESeq2 object from phyloseq object (common taxa)
  rootYxS<-phyloseq_to_deseq2(rootphylo,~1) # make DESeq2 object from phyloseq object (common taxa)
  leafYxS$YxS <- factor(paste0(leafYxS$Harvested,".",leafYxS$Site)) # Make separate "grouping factor" with Year and Site
  rootYxS$YxS <- factor(paste0(rootYxS$Harvested,".",rootYxS$Site)) # Make separate "grouping factor" with Year and Site
  design(leafYxS) <- ~YxS+Genotype+Age # define experimental design: interested in Year x Site averaged over genotypes, ages
  design(rootYxS) <- ~YxS+Genotype+Age # define experimental design: interested in Year x Site averaged over genotypes, ages
  leafYxS = estimateSizeFactors(leafYxS, geoMeans = apply(counts(leafYxS), 1, gm_mean)) # Estimate size factors
  rootYxS = estimateSizeFactors(rootYxS, geoMeans = apply(counts(rootYxS), 1, gm_mean)) # Estimate size factors
  leafYxS<-DESeq(leafYxS,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) # Store differential abundance analysis for this level
  rootYxS<-DESeq(rootYxS,test="Wald",fitType="parametric",parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) # Store differential abundance analysis for this level
  save(leafYxS,file=paste0("foldchange/DESeq2_datasets/dds_leaf_",level,"_YxS.RData")) # save to file
  save(rootYxS,file=paste0("foldchange/DESeq2_datasets/dds_root_",level,"_YxS.RData")) # save to file
  ## Get all pairwise Genotype contrasts  at each site
  for (site in SiteList) { 
      grp1<-paste0("2011.",site)
        grp2<-paste0("2012.",site)
        # print(paste0(grp1,".",grp2)) # debugging line to make sure all pairwise comparisons are represented at each site
        # get site-specific log2 fold change between years at this site
        ## roots first:
        results(rootYxS,contrast=c('YxS',grp1,grp2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
          mutate(Taxon=row.names(.),Term='YxS',Contrast=paste0(grp1,"_",grp2),Site=site,Organ='root',Level=level) %>%
          rbind(.,l2fc.int) -> l2fc.int
        ## then leaves:
        results(leafYxS,contrast=c('YxS',grp1,grp2),parallel=TRUE,BPPARAM=MulticoreParam(workers=5)) %>% as.data.frame %>% 
          mutate(Taxon=row.names(.),Term='YxS',Contrast=paste0(grp1,"_",grp2),Site=site,Organ='leaf',Level=level) %>%
          rbind(.,l2fc.int) -> l2fc.int
  }
}
rm(site,level,grp1,grp2,leafYxS,rootYxS,leafphylo,rootphylo)

####### Add site-specific (interaction term) log2 fold change results to main effects #######
l2fc<-rbind(l2fc,l2fc.int) %>% 
  mutate(Taxon=factor(Taxon),Term=factor(Term),Contrast=factor(Contrast),Site=factor(Site),Organ=factor(Organ),Level=factor(Level))

save(l2fc,file="foldchange/log2foldchange.RData")
rm(l2fc.int)

####### Adjust p-values again for multiple comparisons #######
# For each OTU, there is a p-value for each contrast between each pair of Genotypes, Sites, Ages, etc. 
# P-value for each of these contrasts has already been corrected for multiple comparisons (i.e., 1 test for each OTU)
## Now, get "effect-level" p-values by adjusting again for multiple comparisons within each Term (i.e., 1 test for each pairwise group contrast within each Term)
### This correction asks whether all of the log2 fold changes between group pairs are simultaneously nonzero 
group_by(l2fc,Taxon,Term,Organ) %>% mutate(padjadj=p.adjust(padj,method='BH')) %>%
  ungroup %>% as.data.frame -> l2fc

## Make Level an ordered factor
l2fc$Level<- ordered(l2fc$Level, levels = c("otu", "fam", "ord","cla","phy"))
summary(l2fc$Level)

save(l2fc,file="foldchange/log2foldchange.RData") # save again

####### Median effect sizes on OTU abundances: fold changes  #######

filter(l2fc,padjadj<2) %>% # only look at significant enrichments/depletions
  group_by(Level,Organ,Term) %>% 
  # translate median log2 fold change into fold change; use abs() because we are interested in effect magnitude not direction
  summarize(medianFC=2^median(abs(log2FoldChange),na.rm=TRUE)) %>% 
  ungroup %>% as.data.frame -> medianFCs

filter(medianFCs,Level=='otu')
"   Level Organ     Term  medianFC
1    otu  leaf      Age  2.640838
2    otu  leaf      AxS  7.812085
3    otu  leaf Genotype  3.186153
4    otu  leaf      GxS 15.518764
5    otu  leaf     Site  5.918921
6    otu  leaf     Year  2.272452
7    otu  leaf      YxS  3.688281
8    otu  root      Age  2.226692
9    otu  root      AxS  3.285965
10   otu  root Genotype  2.148515
11   otu  root      GxS  7.300904
12   otu  root     Site  2.994312
13   otu  root     Year  1.843555
14   otu  root      YxS  2.387859 "

####### Fig. 3b: Boxplot: SIGNIFICANT effect sizes (fold changes for each NBM term at all levels) #######
pdf(file="plots/Fig_3b_boxplot_effect_sizes.pdf",width=19,height=12)
update_geom_defaults("point", list(colour = NULL))
filter(l2fc,padjadj<0.05) %>%
  mutate(Term2=paste0("\n",Term,"\n")) %>%
  mutate(Term2=ifelse(Term=='Genotype','\nGeno.\n',Term2)) %>%
  arrange(Level) %>%
  ggplot(.,aes(x=Level,y=abs(log2FoldChange),color=Organ))+
  geom_boxplot(size=1.5,width=1,outlier.colour=NULL,position=position_dodge(w=1),outlier.size=3)+
  facet_grid(.~Term2)+
  coord_cartesian(ylim=c(0, 11))+
  ylab("Fold changes attributed to\neach source of variation\n")+
  scale_color_manual(values=c("forest green","grey"))+
  scale_y_continuous(labels=function(y) format(2^y,digits=2),breaks=c(0,1,2,4,6,8,10))+ # transform from log2 fold change to plain fold change
  scale_x_discrete(breaks=c("otu","fam","ord","cla","phy"),
                   labels=c("OTU","Family","Order","Class","Phylum"))+
  theme_classic()+
  theme(panel.margin = unit(1.35, "lines"))+
  theme(strip.text.x=element_text(size=34,face="bold",lineheight=0.4),strip.background=element_rect(fill="gray95",color="gray95"))+
  theme(axis.text.x = element_text(size=28,angle=55,vjust=1,hjust=1,lineheight=0.2,face="bold"))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=36,face="bold"))+
  theme(legend.title= element_blank(),legend.text=element_text(size=36,face="bold"))+
  theme(axis.text.y=element_text(size=28,face="bold"))+
  theme(legend.background = element_rect(fill="gray95", size=.8))+
  theme(legend.key.height=unit(2.5,"lines"),legend.key.width=unit(2.5,"lines"),legend.position=c(0.9,0.8))
update_geom_defaults("point", list(colour = "black"))
dev.off()

####### Fig. 5b: Dotplot: OTU fold changes due to genotype #######
load('intermediate_data/OTU_tax_table.RData') # load taxonomic info

pdf(file="plots/Fig_5b_Dotplot_OTU_genotype_foldchange_by_Phylum.pdf",height=10,width=15)
set.seed(777) # for reproducible jitter
filter(l2fc,Level=='otu' & Term=='Genotype' & padjadj<0.05) %>%
  merge(.,OTUtax,by.x="Taxon",by.y="OTU_ID") %>% 
  ggplot(.,aes(x=Phylum,y=abs(log2FoldChange),color=Organ,fill=Organ))+
  geom_point(size=5,position=position_jitterdodge(jitter.width=0.7,dodge.width=0.6),alpha=0.4)+
  scale_color_manual(values=c("forest green","dark grey"),guide=FALSE)+
  scale_fill_manual(values=c("forest green","dark grey"),guide=FALSE)+
  ylab("Fold change attributed\nto host genotype\n")+
  xlab("Phylum")+
  scale_y_continuous(labels=function(y) format(2^y,digits=2),breaks=c(0,1,2,3))+
  coord_cartesian(ylim=c(0,3.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size=34,face="bold",angle=65,hjust=0.8,vjust=0.85))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=32,face="bold"))+
  theme(legend.title= element_blank(),legend.text=element_text(size=36,face="bold"))+
  theme(axis.text.y=element_text(size=32,face="bold"),legend.position=c(0.5,0.9))+
  theme(legend.background = element_rect(fill="gray95", size=.8))+
  theme(legend.key.height=unit(2,"line"),legend.key.width=unit(2,"line"))
dev.off()

####### Supplementary Dataset 1 #######

filter(l2fc,padjadj<0.05) %>% # only keep significant enrichment/depletions
  select(Taxon,Organ,Level,Term) %>% 
  unique() %>% # each taxon represented no more than once per term
  mutate(TaxLev=paste0(Taxon,Level)) %>% # to merge with relative abundance data
  select(-Taxon,-Level) %>%
  merge(.,RelAbund,by='TaxLev') %>%
  mutate(RelAbund_5sites=ifelse(Organ=='leaf',RA_leaf5,RA_root5), # give more descriptive names and organ-specific data:
         RelAbund_3sites=ifelse(Organ=='leaf',RA_leaf3,RA_root3),
         RelAbund_Jam=ifelse(Organ=='leaf',RA_leaf_Jam,RA_root_Jam),
         RelAbund_Mah=ifelse(Organ=='leaf',RA_leaf_Mah,RA_root_Mah),
         RelAbund_Mil=ifelse(Organ=='leaf',RA_leaf_Mil,RA_root_Mil),
         RelAbund_Par=ifelse(Organ=='leaf',RA_leaf_Par,RA_root_Par),
         RelAbund_Sil=ifelse(Organ=='leaf',RA_leaf_Sil,RA_root_Sil),
         RelAbund_Age2=ifelse(Organ=='leaf',RA_leaf3_age2,RA_root3_age2),
         RelAbund_Age3=ifelse(Organ=='leaf',RA_leaf3_age3,RA_root3_age3),
         RelAbund_Age4=ifelse(Organ=='leaf',RA_leaf3_age4,RA_root3_age4),
         RelAbund_endog=ifelse(Organ=='leaf',RA_leaf3_endog,RA_root3_endog)) %>%
  mutate(Level=plyr::mapvalues(Level,from=c('otu','fam','ord','cla','phy'),
                               to=c('OTU','Family','Order','Class','Phylum'))) %>%
  select(-starts_with('RA_'),-TaxLev) %>% # get rid of unused columns
  arrange(Term,Level,Organ) %>%
  write.table(.,file="Supplementary_Dataset_1_all_NBM_hits.txt",sep='\t',col.names=TRUE,row.names=FALSE)

####### Load taxon relative abundance data and merge with LRT results #######
for (level in c('otu','fam','ord','cla','phy')){
  load(paste0("rel_abund/RelAbund_withEndog_",level,".RData"))
}

# combine into single dataframe:
RelAbund<-rbind(RelAbund.withEndog.otu,RelAbund.withEndog.fam[,1:31],RelAbund.withEndog.ord,RelAbund.withEndog.cla,RelAbund.withEndog.phy)
rm(RelAbund.withEndog.otu,RelAbund.withEndog.fam,RelAbund.withEndog.ord,RelAbund.withEndog.cla,RelAbund.withEndog.phy)
RelAbund$TaxLev<-paste0(RelAbund$Taxon,RelAbund$Level) # for merging purposes

filter(lrt.all,padj<0.05) %>% # significant LRTs only
  mutate(TaxLev=paste0(Taxon,Level)) %>%
  merge(.,select(RelAbund,TaxLev,RA_root3,RA_leaf3),by='TaxLev') %>%
  mutate(RA3=ifelse(Organ=='root',RA_root3,RA_leaf3)) %>% # give rel. abundances in leaves and roots the same name
  select(-starts_with('RA_'),-TaxLev) %>%
  group_by(Organ,Level,Term) %>% # look separately at each Term for each Level for leaves and roots
  summarize(RelAbund=sum(RA3)) %>%  # add up relative abundances of all taxa with significant LRT
  ungroup %>% as.data.frame -> weighted.NBM.LRT
  
####### Figure 3a: Barplot: WEIGHTED proportions of taxa predicted by NBM LRTs #######

pdf(file="plots/Fig_3a_Barplot_weighted_proportions_NBM_LRT.pdf",width=19,height=12)
mutate(weighted.NBM.LRT,Term2=paste0('\n',Term,'\n')) %>% 
  mutate(Level=ordered(Level,levels=c('otu','fam','ord','cla','phy'))) %>%
  arrange(Level) %>%
  ggplot(.,aes(x=Level,y=RelAbund,fill=Organ))+
  geom_bar(stat="identity",position="dodge",width=1)+
  facet_grid(.~Term2)+
  ylab("Weighted percentage of taxa responsive\nto each source of variation\n")+
  scale_fill_manual(values=c("forest green","grey"))+
  scale_x_discrete(breaks=c("otu","fam","ord","cla","phy"),
                   labels=c("OTU","Family","Order","Class","Phylum"))+
  theme_classic()+scale_y_continuous(labels=function(y) 100*y)+
  theme(panel.margin = unit(1.35, "lines"))+
  theme(strip.text.x=element_text(size=34,face="bold",lineheight=0.3),strip.background=element_rect(fill="gray95",color="gray95"))+
  theme(axis.text.x = element_text(size=28,angle=55,vjust=1,hjust=1,lineheight=0.2,face="bold"))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=36,face="bold"))+
  theme(legend.title= element_blank(),legend.text=element_text(size=36,face="bold"))+
  theme(axis.text.y=element_text(size=28,face="bold"),legend.position=c(0.4,0.85))+
  theme(legend.background = element_rect(fill="gray95", size=.8))+
  theme(legend.key.height=unit(2.5,"lines"),legend.key.width=unit(2.5,"lines"))
dev.off()

####### Supplementary Dataset 1 #######

filter(lrt.all,padj<0.05) %>% # only keep significant LRT hits
  select(Taxon,Organ,Level,Term,padj) %>% 
  mutate(TaxLev=paste0(Taxon,Level)) %>% # to merge with relative abundance data
  select(-Taxon,-Level) %>%
  merge(.,RelAbund,by='TaxLev') %>%
  mutate(RelAbund_5sites=ifelse(Organ=='leaf',RA_leaf5,RA_root5), # give more descriptive names and organ-specific data:
         RelAbund_3sites=ifelse(Organ=='leaf',RA_leaf3,RA_root3),
         RelAbund_Jam=ifelse(Organ=='leaf',RA_leaf_Jam,RA_root_Jam),
         RelAbund_Mah=ifelse(Organ=='leaf',RA_leaf_Mah,RA_root_Mah),
         RelAbund_Mil=ifelse(Organ=='leaf',RA_leaf_Mil,RA_root_Mil),
         RelAbund_Par=ifelse(Organ=='leaf',RA_leaf_Par,RA_root_Par),
         RelAbund_Sil=ifelse(Organ=='leaf',RA_leaf_Sil,RA_root_Sil),
         RelAbund_Age2=ifelse(Organ=='leaf',RA_leaf3_age2,RA_root3_age2),
         RelAbund_Age3=ifelse(Organ=='leaf',RA_leaf3_age3,RA_root3_age3),
         RelAbund_Age4=ifelse(Organ=='leaf',RA_leaf3_age4,RA_root3_age4),
         RelAbund_endog=ifelse(Organ=='leaf',RA_leaf3_endog,RA_root3_endog)) %>%
  mutate(Level=plyr::mapvalues(Level,from=c('otu','fam','ord','cla','phy'),
                               to=c('OTU','Family','Order','Class','Phylum'))) %>%
  select(-starts_with('RA_'),-TaxLev) %>% # get rid of unused columns
  arrange(Term,Level,Organ) %>%
  write.table(.,file="Supplementary_Dataset_1_all_LRT_hits.txt",sep='\t',col.names=TRUE,row.names=FALSE)

####### Figure 6d: Barplot: weighted proportion of OTUs differing among genotypes AT EACH SITE ####### 

pdf(file="plots/Fig_6d_barplot_OTU_genotypeXsite.pdf")
filter(l2fc,padjadj<0.05,Level=='otu',Term=='GxS') %>% # significant OTU enrichments/depletions only
  select(Taxon,Term,Site,Organ,RA3) %>%
  unique() %>% # make sure each OTU is not represented >1 time at each site 
  group_by(Organ,Site) %>%
  summarize(siteRA=sum(RA3)) %>% # sum relative abundance in each site
  ungroup %>% as.data.frame %>%
  ggplot(.,aes(x=Organ,color=Site,fill=Site,y=siteRA))+
  geom_bar(stat='identity',position='dodge')+
  scale_color_manual(values=sitePalette)+
  scale_fill_manual(values=sitePalette)+
  ylab("Weighted percentage of\nOTUs differentially abundant\nbetween genotypes\n")+
  theme_classic()+scale_y_continuous(labels=function(y) 100*y)+
  theme(axis.text.x = element_text(size=36,face="bold",color=c('forest green','dark grey')))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=38,face="bold"))+
  theme(legend.title= element_text(size=34,face="bold"),legend.text=element_text(size=28,face="bold"))+
  theme(axis.text.y=element_text(size=30,face="bold"),legend.position=c(0.8,0.8))+
  theme(legend.background = element_rect(fill="gray95", size=.8))+
  theme(legend.key.height=unit(2.5,"lines"),legend.key.width=unit(2.5,"lines"))
dev.off()

####### Figure 4d: Plot Phylum changes, Site*Age, leaf #######

filter(lrt.all,Organ=='leaf',Term=='AxS',Level=='phy',padj<0.05) # Actinobacteria, Armatimonadetes, and Verrucomicrobia have significant Age*Site interaction

pdf(file="plots/Fig_4d_Phyla_AxS_leaf.pdf",width=9,height=11)
filter(l2fc,Organ=='leaf',Term=='AxS',Level=='phy',Taxon%in%c('Actinobacteria','Armatimonadetes','Verrucomicrobia')) %>%
  mutate(Contrast=factor(gsub("[^0-9]+","",Contrast))) %>%  # simplify Contrast info
  mutate(Contrast=paste0('ages',Contrast))%>%
  rename(Age2=baseMean) %>%
  select(-stat,-pvalue,-padj,-padjadj,-lfcSE) %>% # only keep relevant columns
  group_by(Taxon,Site) %>%
  spread(log2FoldChange,key=Contrast) %>%
  mutate(Age3=Age2*(2^ages23), Age4=Age3*(2^ages34)) %>% # calculate mean at each age from starting mean and log2 fold changes
  gather(key=Age,value=Abundance,Age2,Age3,Age4) %>% mutate(Age=gsub("[^0-9]","",Age)) %>%
  ungroup %>% as.data.frame %>% 
  mutate(Taxon=plyr::mapvalues(Taxon,from=c('Actinobacteria','Armatimonadetes','Verrucomicrobia'),
                               to=c('Actino.','Armatim.','Verruco.'))) %>%
  ggplot(.,aes(x=Age,y=Abundance,color=Site,group=Site))+
  geom_point(size=5,alpha=1,position=position_dodge(w=0.2))+
  geom_line(size=2,position=position_dodge(w=0.2))+
  facet_grid(Taxon~.,scales="free_y")+
  scale_color_manual(values=sitePalette)+
  ylab("Mean abundance in leaves")+xlab("Plant age (years)")+
  theme_minimal()+theme(legend.background=element_rect(fill="gray90",color="gray90"))+
  theme(strip.text.y = element_text(size=32,face="bold"))+theme(legend.position="top",legend.key.height=unit(1.5,"lines"),legend.key.width=unit(1.5,"lines"))+
  theme(axis.title.y = element_text(size=36,face="bold"),axis.text.y=element_text(size=22,face="bold"))+
  theme(axis.title.x = element_text(size=36,face="bold"),axis.text.x=element_text(size=30,face="bold"))+
  theme(legend.title = element_text(size=30,face="bold"),legend.text=element_text(size=28,face="bold"))+
  theme(panel.margin = unit(1.5, "lines"),strip.background=element_rect(fill="gray90",color="gray90"))
dev.off()

####### FIGURE S9b: Re-create Fig. 3b with these 29 OTUs to show that they are representative of the rest of the OTUs #######
load('succession/l2fc_nonEndog_results.RData') # load fold change results for the "non-endogenous" OTUs

l2fc.otu<-filter(l2fc,Level=='otu') %>% mutate('Origin'='Natural') %>% select(-RA3,-Level)
mutate(l2fc.nonEndog,Origin='Unnatural') %>% rename(Taxon=OTU_ID) %>%
  rbind(.,l2fc.otu) %>% filter(padjadj<0.05) -> l2fc.otu

pdf(file="plots/Fig_S9b_foldchange_unnatural.pdf",width=11,height=9)
ggplot(l2fc.otu,aes(x=Term,y=abs(log2FoldChange),color=Origin,fill=Origin))+
  facet_grid(.~Organ)+
  geom_point(size=3,position=position_jitterdodge(),alpha=0.3)+
  scale_y_continuous(labels=function(y) format(2^y,digits=2),breaks=c(0,1,2,4,6,8,10,12,16))+ # transform from log2 fold change to plain fold change
  theme_classic()+
  ylab("Fold changes attributed to\neach source of variation")+
  scale_color_manual(values=c('dark grey','royal blue'))+
  scale_fill_manual(values=c('dark grey','royal blue'))+
  theme(panel.margin = unit(1.35, "lines"))+
  theme(strip.text.x=element_text(size=34,face="bold",lineheight=0.4),strip.background=element_rect(fill="gray95",color="gray95"))+
  theme(axis.text.x = element_text(size=24,angle=55,vjust=1,hjust=1,lineheight=0.2,face="bold"))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=36,face="bold"))+
  theme(legend.title= element_blank(),legend.text=element_text(size=30,face="bold"))+
  theme(axis.text.y=element_text(size=24))+
  theme(legend.background = element_rect(fill="gray95", size=.8),legend.position='top')+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(1.5,"lines"))
dev.off()


####### Save image #######
save.image(paste0("foldchange/image_",date(),".RData"))



