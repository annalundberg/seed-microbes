### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
#####  Analyses to explore the influence of transplant from greenhouse on microbiome change with plant age
####### Clear workspace ########
rm(list=ls())
####### Load source file #######
source("ecotypes_source.R")

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
[5] RSQLite_1.0.0        chron_2.3-47         RColorBrewer_1.1-2   minqa_1.2.4         
[9] colorspace_1.2-6     XML_3.98-1.4         zlibbioc_1.16.0      xtable_1.8-2        
[13] annotate_1.48.0      mgcv_1.8-12          nnet_7.3-12          survival_2.38-3     
[17] RJSONIO_1.3-0        magrittr_1.5         nlme_3.1-126         MASS_7.3-45         
[21] foreign_0.8-66       tools_3.2.3          data.table_1.9.6     stringr_1.0.0       
[25] munsell_0.4.3        locfit_1.5-9.1       cluster_2.0.3        AnnotationDbi_1.32.3
[29] lambda.r_1.1.7       ade4_1.7-4           nloptr_1.0.4         biom_0.3.12         
[33] igraph_1.0.1         gtable_0.2.0         codetools_0.2-14     multtest_2.26.0     
[37] DBI_0.3.1            R6_2.1.2             gridExtra_2.2.1      Hmisc_3.17-2        
[41] futile.options_1.0.0 stringi_1.0-1        geneplotter_1.48.0   rpart_4.1-10        
[45] acepack_1.3-3.3 "

####### Load intermediate data files: Phyloseq object with FIELD samples & with ALL samples #######
load("intermediate_data/fieldEco_cleaned_CNC.RData")
load("intermediate_data/fullEco_cleaned_CNC.RData")

# rename to something shorter for convenience:
fieldEco<-fieldEco.nobadOTUs.highcoverage.thresholded.CNC; rm(fieldEco.nobadOTUs.highcoverage.thresholded.CNC)
fullEco<-fullEco.nobadOTUs.highcoverage.thresholded.CNC; rm(fullEco.nobadOTUs.highcoverage.thresholded.CNC)

####### List of tests #######
# 1. How many OTUs are present in experimental but not endogenous plants, and how common are they?
# 2. Does dissimilarity among sites increase with age in roots?
# 2a. in leaves?
# 3. Do roots become more similar to endogenous roots with age?
# 3a. in leaves?

####### How many OTUs are present in experimental but not endogenous plants? #######
## Make subsets of endogenous/experimental  samples
EndogAndSoilOnly<-subset_samples(fieldEco,Genotype%in%c('endog','soil')) # separate out Endogenous plants
EndogAndSoilOnly<-prune_taxa(taxa_sums(EndogAndSoilOnly)>0,EndogAndSoilOnly) # remove OTUs with no observations
ntaxa(EndogAndSoilOnly) # 3689 OTUs remain
nonEndog.OTUs<-(setdiff(taxa_names(fieldEco),taxa_names(EndogAndSoilOnly)))
length(nonEndog.OTUs)
# 29 OTUs not seen in endogenous samples
save(nonEndog.OTUs,file="succession/non_endogenous_OTUs.RData")

####### These "unnatural" OTUs: were they evenly distributed in the experiment? #######
nonEndog.root.phyloseq<-prune_taxa(nonEndog.OTUs,fieldEco) %>% # only keep these 29 OTUs
  subset_samples(Type=='root' & Age!='endog' & Site%in%c('Jam','Mah','Sil')) %>% # only keep root samples from experimental plants at 3 main sites
  prune_samples(sample_sums(.)>0,.) # throw out any samples that don't have any of these 29 OTUs
nsamples(nonEndog.root.phyloseq)  # 307 of 310 root samples left

####### Create DESeq2 datasets for MAIN EFFECTS (no interaction terms) #######

## See DESeq2 vignette for explanation of why including interaction terms alters interpretation of main effects (justification for doing 2 separate models: main effects only and interactions only)
dds.nonEndog.root.main<-phyloseq_to_deseq2(nonEndog.root.phyloseq,~Site+Genotype+Age+Harvested)
dds.nonEndog.root.main = estimateSizeFactors(dds.nonEndog.root.main, geoMeans = apply(counts(dds.nonEndog.root.main), 1, gm_mean))
dds.nonEndog.root.main<-DESeq(dds.nonEndog.root.main,test="Wald",fitType="parametric")

# Do the same for leaves:
nonEndog.leaf.phyloseq<-prune_taxa(nonEndog.OTUs,fieldEco) %>% # only keep these 29 OTUs
  subset_samples(Type=='leaf' & Age!='endog' & Site%in%c('Jam','Mah','Sil')) %>% # only keep root samples from experimental plants at 3 main sites
  prune_samples(sample_sums(.)>0,.) # throw out any samples that don't have any of these 29 OTUs
nsamples(nonEndog.leaf.phyloseq)  # 197 of 306 leaf samples left

# Create DESeq2 object
dds.nonEndog.leaf.main<-phyloseq_to_deseq2(nonEndog.leaf.phyloseq,~Site+Genotype+Age+Harvested)
dds.nonEndog.leaf.main = estimateSizeFactors(dds.nonEndog.leaf.main, geoMeans = apply(counts(dds.nonEndog.leaf.main), 1, gm_mean))
dds.nonEndog.leaf.main<-DESeq(dds.nonEndog.leaf.main,test="Wald",fitType="parametric")

###### otu: save main-effects DESeq objects #######
save(dds.nonEndog.leaf.main, file="succession/dds_nonEndog_leaf_main.RData")
save(dds.nonEndog.root.main,file="succession/dds_nonEndog_root_main.RData")

####### otu: likelihood ratio tests for main effects #######
# LRT for Genotype: 
nbinomLRT(dds.nonEndog.root.main,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root') -> lrt.nonEndog
nbinomLRT(dds.nonEndog.leaf.main,reduced=~Site+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog

# LRT for Site: 
nbinomLRT(dds.nonEndog.root.main,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog
nbinomLRT(dds.nonEndog.leaf.main,reduced=~Genotype+Age+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog

# LRT for Age: 
nbinomLRT(dds.nonEndog.root.main,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog
nbinomLRT(dds.nonEndog.leaf.main,reduced=~Site+Genotype+Harvested) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog

# LRT for Year: 
nbinomLRT(dds.nonEndog.root.main,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='root')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog
nbinomLRT(dds.nonEndog.leaf.main,reduced=~Site+Genotype+Age) %>% results %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Year',Organ='leaf')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog
####### otu: likelihood ratio tests for interaction effects #######
## LRT for Genotype:Site
nbinomLRT(dds.nonEndog.root.main,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.nonEndog.root.main)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='root')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog
nbinomLRT(dds.nonEndog.leaf.main,full=~Site+Genotype+Age+Harvested+Genotype:Site,reduced=design(dds.nonEndog.leaf.main)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='GxS',Organ='leaf')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog

## LRT for Age:Site
nbinomLRT(dds.nonEndog.root.main,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.nonEndog.root.main)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='root')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog
nbinomLRT(dds.nonEndog.leaf.main,full=~Site+Genotype+Age+Harvested+Age:Site,reduced=design(dds.nonEndog.leaf.main)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='AxS',Organ='leaf')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog

## LRT for Year:Site
nbinomLRT(dds.nonEndog.root.main,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.nonEndog.root.main)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='root')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog
nbinomLRT(dds.nonEndog.leaf.main,full=~Site+Genotype+Age+Harvested+Harvested:Site,reduced=design(dds.nonEndog.leaf.main)) %>%
  results %>% as.data.frame %>% mutate(OTU_ID=row.names(.),Term='YxS',Organ='leaf')  %>% rbind(.,lrt.nonEndog) -> lrt.nonEndog
####### save non-endogenous OTU likelihood ratio test results #######
save(lrt.nonEndog,file="succession/LRT_nonEndog_results.RData")

####### otu: log2 fold changes due to main effects #######
## Get all pairwise Genotype contrasts for roots:
results(dds.nonEndog.root.main,contrast=c("Genotype","JAM","MAH")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='JAM_MAH') -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","JAM","MIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='JAM_MIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","JAM","PAR")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='JAM_PAR') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","JAM","SIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='JAM_SIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","MAH","MIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='MAH_MIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","MAH","PAR")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='MAH_PAR') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","MAH","SIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='MAH_SIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","MIL","PAR")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='MIL_PAR') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","MIL","SIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='MIL_SIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Genotype","PAR","SIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='root',Contrast='PAR_SIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
# Get all pairwise Genotype contrasts for leaves:
results(dds.nonEndog.leaf.main,contrast=c("Genotype","JAM","MAH")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='JAM_MAH') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","JAM","MIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='JAM_MIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","JAM","PAR")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='JAM_PAR') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","JAM","SIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='JAM_SIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","MAH","MIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='MAH_MIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","MAH","PAR")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='MAH_PAR') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","MAH","SIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='MAH_SIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","MIL","PAR")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='MIL_PAR') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","MIL","SIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='MIL_SIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Genotype","PAR","SIL")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Genotype',Organ='leaf',Contrast='PAR_SIL') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog

# Get all pairwise Site contrasts for roots:
results(dds.nonEndog.root.main,contrast=c("Site","Jam","Mah")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root',Contrast='Jam_Mah') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Site","Jam","Sil")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root',Contrast='Jam_Sil') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Site","Mah","Sil")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='root',Contrast='Mah_Sil') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
# Get all pairwise Site contrasts for leaves:
results(dds.nonEndog.leaf.main,contrast=c("Site","Jam","Mah")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf',Contrast='Jam_Mah') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Site","Jam","Sil")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf',Contrast='Jam_Sil') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Site","Mah","Sil")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Site',Organ='leaf',Contrast='Mah_Sil') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog

# Get all pairwise Age contrasts for roots:
results(dds.nonEndog.root.main,contrast=c("Age","2","3")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root',Contrast='2_3') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Age","2","4")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root',Contrast='2_4') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.main,contrast=c("Age","3","4")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='root',Contrast='3_4') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
# Get all pairwise Age contrasts for leaves:
results(dds.nonEndog.leaf.main,contrast=c("Age","2","3")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf',Contrast='2_3') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Age","2","4")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf',Contrast='2_4') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Age","3","4")) %>% 
  as.data.frame %>% mutate(OTU_ID=row.names(.),Term='Age',Organ='leaf',Contrast='3_4') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog

# Get Year effect for roots and leaves:
results(dds.nonEndog.root.main,contrast=c("Harvested","2011","2012")) %>% as.data.frame %>% 
  mutate(OTU_ID=row.names(.),Term='Year',Organ='root',Contrast='2011_2012') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.main,contrast=c("Harvested","2011","2012")) %>% as.data.frame %>% 
  mutate(OTU_ID=row.names(.),Term='Year',Organ='leaf',Contrast='2011_2012') %>% rbind(.,l2fc.nonEndog) -> l2fc.nonEndog

## Add column to note that the main effects apply averaged across all sites (to differentiate from interaction term)
l2fc.nonEndog$Site<-'all'

####### otu: site-specific log2 fold changes (interaction terms) #######
####### log2 fold changes: Genotype x Site #######
## roots: make DESeq2 object to look at GxS while controlling for Age, Year
dds.nonEndog.root.GxS<-phyloseq_to_deseq2(nonEndog.root.phyloseq,~1) # make DESeq2 object
dds.nonEndog.root.GxS$GxS <- factor(paste0(dds.nonEndog.root.GxS$Genotype,".",dds.nonEndog.root.GxS$Site)) # Make separate "grouping factor" with Genotype and Site
design(dds.nonEndog.root.GxS) <- ~GxS+Age+Harvested # define experimental design: interested in Geno x Site averaged over ages, years
dds.nonEndog.root.GxS = estimateSizeFactors(dds.nonEndog.root.GxS, geoMeans = apply(counts(dds.nonEndog.root.GxS), 1, gm_mean))
dds.nonEndog.root.GxS<-DESeq(dds.nonEndog.root.GxS,test="Wald",fitType="parametric")
## leaves: make DESeq2 object to look at GxS while controlling for Age, Year
dds.nonEndog.leaf.GxS<-phyloseq_to_deseq2(nonEndog.leaf.phyloseq,~1) # make DESeq2 object
dds.nonEndog.leaf.GxS$GxS <- factor(paste0(dds.nonEndog.leaf.GxS$Genotype,".",dds.nonEndog.leaf.GxS$Site)) # Make separate "grouping factor" with Genotype and Site
design(dds.nonEndog.leaf.GxS) <- ~GxS+Age+Harvested # define experimental design: interested in Geno x Site averaged over ages, years
dds.nonEndog.leaf.GxS = estimateSizeFactors(dds.nonEndog.leaf.GxS, geoMeans = apply(counts(dds.nonEndog.leaf.GxS), 1, gm_mean))
dds.nonEndog.leaf.GxS<-DESeq(dds.nonEndog.leaf.GxS,test="Wald",fitType="parametric")

## Get all pairwise Genotype contrasts for roots at each site
GenoList<-c('JAM','MAH','MIL','PAR','SIL')
SiteList<-c('Jam','Mah','Sil')
for (i in 1:length(SiteList)) { # i = index of site
  site<-SiteList[i] 
  for (j in 1:(length(GenoList)-1)) { # j = index of 1st genotype
    geno1<-paste0(GenoList[j],".",site)
    for (k in (j+1):length(GenoList)) { # k = index of 2nd genotype
      geno2<-paste0(GenoList[k],".",site)
      # print(paste0(geno1,".",geno2)) # debugging line to make sure all pairwise comparisons are represented at each site
      # get site-specific log2 fold change between this pair of genotypes at site i
      ## roots first:
      results(dds.nonEndog.root.GxS,contrast=c('GxS',geno1,geno2)) %>% as.data.frame %>% 
        mutate(OTU_ID=row.names(.),Term='GxS',Contrast=paste0(geno1,"_",geno2),Site=site,Organ='root') %>%
        rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
      ## then leaves:
      results(dds.nonEndog.leaf.GxS,contrast=c('GxS',geno1,geno2)) %>% as.data.frame %>% 
        mutate(OTU_ID=row.names(.),Term='GxS',Contrast=paste0(geno1,"_",geno2),Site=site,Organ='leaf') %>%
        rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
    }
  }
}
rm(site,geno1,geno2,i,j,k) # clean up

####### log2 fold changes: Age x Site #######
## roots: make DESeq2 object to look at AxS while controlling for Genotype, Year
dds.nonEndog.root.AxS<-phyloseq_to_deseq2(nonEndog.root.phyloseq,~1) # make DESeq2 object
dds.nonEndog.root.AxS$AxS <- factor(paste0(dds.nonEndog.root.AxS$Age,".",dds.nonEndog.root.AxS$Site)) # Make separate "grouping factor" with Age and Site
design(dds.nonEndog.root.AxS) <- ~AxS+Genotype+Harvested # define experimental design: interested in Age x Site averaged over genotypes, years
dds.nonEndog.root.AxS = estimateSizeFactors(dds.nonEndog.root.AxS, geoMeans = apply(counts(dds.nonEndog.root.AxS), 1, gm_mean))
dds.nonEndog.root.AxS<-DESeq(dds.nonEndog.root.AxS,test="Wald",fitType="parametric")
## leaves: make DESeq2 object to look at AxS while controlling for Genotype, Year
dds.nonEndog.leaf.AxS<-phyloseq_to_deseq2(nonEndog.leaf.phyloseq,~1) # make DESeq2 object
dds.nonEndog.leaf.AxS$AxS <- factor(paste0(dds.nonEndog.leaf.AxS$Age,".",dds.nonEndog.leaf.AxS$Site)) # Make separate "grouping factor" with Age and Site
design(dds.nonEndog.leaf.AxS) <- ~AxS+Genotype+Harvested # define experimental design: interested in Age x Site averaged over genotypes, years
dds.nonEndog.leaf.AxS = estimateSizeFactors(dds.nonEndog.leaf.AxS, geoMeans = apply(counts(dds.nonEndog.leaf.AxS), 1, gm_mean))
dds.nonEndog.leaf.AxS<-DESeq(dds.nonEndog.leaf.AxS,test="Wald",fitType="parametric")

## Get all pairwise Age contrasts for roots at each site
AgeList<-c('2','3','4')
SiteList<-c('Jam','Mah','Sil')
for (i in 1:length(SiteList)) { # i = index of site
  site<-SiteList[i] 
  for (j in 1:(length(AgeList)-1)) { # j = index of 1st Age
    Age1<-paste0(AgeList[j],".",site)
    for (k in (j+1):length(AgeList)) { # k = index of 2nd Age
      Age2<-paste0(AgeList[k],".",site)
      # print(paste0(Age1,".",Age2)) # debugging line to make sure all pairwise comparisons are represented at each site
      # get site-specific log2 fold change between this pair of Agetypes at site i
      ## roots first:
      results(dds.nonEndog.root.AxS,contrast=c('AxS',Age1,Age2)) %>% as.data.frame %>% 
        mutate(OTU_ID=row.names(.),Term='AxS',Contrast=paste0(Age1,"_",Age2),Site=site,Organ='root') %>%
        rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
      ## then leaves:
      results(dds.nonEndog.leaf.AxS,contrast=c('AxS',Age1,Age2)) %>% as.data.frame %>% 
        mutate(OTU_ID=row.names(.),Term='AxS',Contrast=paste0(Age1,"_",Age2),Site=site,Organ='leaf') %>%
        rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
    }
  }
} 
rm(site,Age1,Age2,i,j,k) # clean up

####### log2 fold changes: Year x Site #######
## roots: make DESeq2 object to look at YxS while controlling for Genotype, Age
dds.nonEndog.root.YxS<-phyloseq_to_deseq2(nonEndog.root.phyloseq,~1) # make DESeq2 object
dds.nonEndog.root.YxS$YxS <- factor(paste0(dds.nonEndog.root.YxS$Harvested,".",dds.nonEndog.root.YxS$Site)) # Make separate "grouping factor" with Year and Site
design(dds.nonEndog.root.YxS) <- ~YxS+Genotype+Age # define experimental design: interested in Year x Site averaged over genotypes, ages
dds.nonEndog.root.YxS = estimateSizeFactors(dds.nonEndog.root.YxS, geoMeans = apply(counts(dds.nonEndog.root.YxS), 1, gm_mean))
dds.nonEndog.root.YxS<-DESeq(dds.nonEndog.root.YxS,test="Wald",fitType="parametric")
## leaves: make DESeq2 object to look at YxS while controlling for Genotype, Age
dds.nonEndog.leaf.YxS<-phyloseq_to_deseq2(nonEndog.leaf.phyloseq,~1) # make DESeq2 object
dds.nonEndog.leaf.YxS$YxS <- factor(paste0(dds.nonEndog.leaf.YxS$Harvested,".",dds.nonEndog.leaf.YxS$Site)) # Make separate "grouping factor" with Year and Site
design(dds.nonEndog.leaf.YxS) <- ~YxS+Genotype+Age # define experimental design: interested in Year x Site averaged over genotypes, ages
dds.nonEndog.leaf.YxS = estimateSizeFactors(dds.nonEndog.leaf.YxS, geoMeans = apply(counts(dds.nonEndog.leaf.YxS), 1, gm_mean))
dds.nonEndog.leaf.YxS<-DESeq(dds.nonEndog.leaf.YxS,test="Wald",fitType="parametric")

## Get Year contrasts at each site:
### Roots first:
results(dds.nonEndog.root.YxS,contrast=c('YxS','2011.Jam','2012.Jam')) %>% as.data.frame %>% 
  mutate(OTU_ID=row.names(.),Term='YxS',Contrast='2011.Jam_2012.Jam',Site='Jam',Organ='root') %>%
  rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.YxS,contrast=c('YxS','2011.Mah','2012.Mah')) %>% as.data.frame %>% 
  mutate(OTU_ID=row.names(.),Term='YxS',Contrast='2011.Mah_2012.Mah',Site='Mah',Organ='root') %>%
  rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.root.YxS,contrast=c('YxS','2011.Sil','2012.Sil')) %>% as.data.frame %>% 
  mutate(OTU_ID=row.names(.),Term='YxS',Contrast='2011.Sil_2012.Sil',Site='Sil',Organ='root') %>%
  rbind(.,l2fc.nonEndog) -> l2fc.nonEndog

### Then leaves:
results(dds.nonEndog.leaf.YxS,contrast=c('YxS','2011.Jam','2012.Jam')) %>% as.data.frame %>% 
  mutate(OTU_ID=row.names(.),Term='YxS',Contrast='2011.Jam_2012.Jam',Site='Jam',Organ='leaf') %>%
  rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.YxS,contrast=c('YxS','2011.Mah','2012.Mah')) %>% as.data.frame %>% 
  mutate(OTU_ID=row.names(.),Term='YxS',Contrast='2011.Mah_2012.Mah',Site='Mah',Organ='leaf') %>%
  rbind(.,l2fc.nonEndog) -> l2fc.nonEndog
results(dds.nonEndog.leaf.YxS,contrast=c('YxS','2011.Sil','2012.Sil')) %>% as.data.frame %>% 
  mutate(OTU_ID=row.names(.),Term='YxS',Contrast='2011.Sil_2012.Sil',Site='Sil',Organ='leaf') %>%
  rbind(.,l2fc.nonEndog) -> l2fc.nonEndog

####### Adjust p-values again for multiple comparisons #######
# For each OTU, there is a p-value for each contrast between each pair of Genotypes, Sites, Ages, etc. 
# P-value for each of these contrasts has already been corrected for multiple comparisons (i.e., 1 test for each OTU)
## Now, get "effect-level" p-values by adjusting again for multiple comparisons within each Term (i.e., 1 test for each pairwise group contrast within each Term)
### This correction asks whether all of the log2 fold changes between group pairs are simultaneously nonzero 
group_by(l2fc.nonEndog,OTU_ID,Term,Organ) %>% mutate(padjadj=p.adjust(padj,method='BH')) %>%
  ungroup %>% as.data.frame -> l2fc.nonEndog

####### Save all log2 fold change results #######
save(l2fc.nonEndog,file='succession/l2fc_nonEndog_results.RData')
## Fig. S9b will be made later using this dataframe (see '../foldchange/foldchange.R')

####### FIGURE S9a: Do these 'non Endogenous' OTUs reliably decrease in abundance with age, as expected? #######
pdf(file="plots/FigS9a_nonEndog_OTU_foldchange_Age.pdf",width=9,height=9)
set.seed(777)
filter(l2fc.nonEndog,Term=='Age') %>%
  mutate(sig=ifelse(padjadj<0.05,'sig','ns')) %>%
  filter(Contrast=='2_4') %>% # only show net changes from 2-4 years of age
  ggplot(.,aes(x=OTU_ID,y=log2FoldChange,color=Organ,fill=Organ))+
  geom_point(size=4,position=position_jitterdodge(jitter.width=0.4,dodge.width=0.6),alpha=0.3)+
  geom_point(data=filter(l2fc.nonEndog,Contrast=='2_4',padjadj<0.05),
             aes(x=OTU_ID,y=log2FoldChange,color=Organ,fill=Organ),
             size=6,position=position_jitterdodge(jitter.width=0.4,dodge.width=0.6))+
  geom_hline(yintercept=0,lty=2)+
  scale_color_manual(values=c("forest green","dark grey"))+
  scale_fill_manual(values=c("forest green","dark grey"))+
  ylab("Fold change in abundance\nbetween plant ages 2 and 4")+
  xlab("non-endogenous OTU")+
  scale_y_continuous(breaks=c(-10,-5,-2,0,2,5,10),labels=c('-1024','-32', '-4', '','4','32','1024'))+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_text(size=32,face="bold"),axis.title.y=element_text(size=32,face="bold"))+
  theme(legend.title= element_blank(),legend.text=element_text(size=36,face="bold"))+
  theme(axis.text.y=element_text(size=28))+
  theme(legend.background = element_rect(fill="gray95", size=.8))+
  theme(legend.key.height=unit(2,"line"),legend.key.width=unit(2,"line"))
dev.off()

####### Supp. Dataset 2: taxonomy of putative 'unnatural' OTUs  #######
load('./rel_abund/RelAbund_withEndog_otu.RData')
data.frame("OTU_ID"=names(refseq(nonEndog.root.phyloseq)), # store OTU names
           "RefSeq"=as.character(refseq(nonEndog.root.phyloseq))) %>% # get reference sequence
  merge(.,tax_table(nonEndog.root.phyloseq),by.x='OTU_ID',by.y='row.names') %>% # merge with taxonomic designation
  merge(.,RelAbund.withEndog.otu,by.x='OTU_ID',by.y='Taxon') %>% # merge with calculated relative abundance
  select(-starts_with('RA_soil'),-ends_with('Par'),-ends_with('Mil'),-Level) %>% # pare down dataset to include main experiment only
write.table(.,file="Supplementary_Dataset_2.txt",sep='\t',row.names=FALSE,col.names=TRUE) # write to file

####### How many OTUs in roots of endogenous plants are also found in potting soil / roots of plants grown in potting soil in the greenhouse? #######
GreenhouseOnly<-subset_samples(fullEco,Treatment=='GRHsoil') %>% prune_taxa(taxa_sums(.)>0,.)
intersect(taxa_names(GreenhouseOnly),taxa_names(EndogAndSoilOnly)) %>% length 
# 3117 OTUs found in both Endogenous roots/wild soils AND in potting soil/roots in potting soil sampled in the greenhouse
3117/ntaxa(EndogAndSoilOnly) # 84.5%

####### Load taxon relative abundance data (if not already loaded) #######
for (i in 1:5){
  level<-c('otu','fam','ord','cla','phy')[i]
  load(paste0("rel_abund/RelAbund_withEndog_",level,".RData"))
}

#######  How abundant are the OTUs that are present in experimental but not endogenous plants? ######
filter(RelAbund.withEndog.otu,Taxon%in%nonEndog.OTUs) %>% select(RA_leaf3) %>% sum
# 0.0021 in leaves at the 3 main sites

filter(RelAbund.withEndog.otu,Taxon%in%nonEndog.OTUs) %>% select(RA_root3) %>% sum
# 0.0011 in roots at the 3 main sites 

####### Alpha diversity: compare between potting-soil and field plants #######
# Use non-transformed OTU tables
pdf(file="plots/Fig_S10b_alphadiv_greenhouse_rootsoil.pdf",width=9,height=9)
update_geom_defaults("point", list(colour = NULL))
subset_samples(fullEco,Type!='leaf' & Treatment%in%c('field','GRHsoil') & !(Site%in%c('Mil','Par'))) %>%
  prune_taxa(taxa_sums(.)>0,.) %>% 
  estimate_richness(.,measures=c('Chao1','Shannon')) %>% as.data.frame() %>%
  merge(.,as(sample_data(fullEco),'data.frame'),by.x='row.names',by.y='SampleID') %>%
  melt(.,measure.vars=c('Chao1','Shannon')) %>% 
  ggplot(.,aes(x=Type,y=value,color=Site))+
    facet_wrap(~variable,scales="free")+
    geom_boxplot()+
    theme_classic()+
    scale_color_manual(values=c('forest green',sitePalette),labels=c('greenhouse','Jam','Mah','Sil'))+
  ylab("Alpha-diversity")+
  theme(strip.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(size=24,face="bold"))+
  theme(axis.title.y=element_text(size=26,face="bold"),axis.text.y=element_text(size=20))+
  theme(legend.title=element_text(size=26,face='bold'),legend.text=element_text(size=22,face='bold'))
update_geom_defaults("point", list(colour = 'black'))
dev.off()


####### PERMUTATION TEST PREDICTIONS & STRATEGY #######
## IF succession is still happening, we expect: 
### 1. Sites get more distinct with age (in roots but not leaves)
### 2. Roots but not leaves get more similar to endogenous plants with age
### 3. Roots get less similar to GRH potting soil

#### Approach:
### 1. Get Weighted UniFrac distances for ENTIRE dataset incl. endog, GRH
### 2. Extract data frame with W. UniFrac PCo1-PCo3 for all samples
### 3. Calculate MEAN PCo1-3 for Potting Soil controls, then remove from dataframe
### 4. Calculate MEAN PCo1-3 for Endogenous plants (Leaf and Root) AT EACH SITE, then remove from df
### 5. Fit standard model to experimental leaf/root samples from 3 major sites; extract Site*Age LS means
#### Test prediction: Sites get more distinct with age (in roots but not leaves)
### 6. Calculate variance of LS means for 3 sites in each age group, fit linear functions of age, save Betas
###    (do for both leaves and roots)
#### Test prediction: Roots but not leaves get more similar to endogenous plants with age
### 7. Calculate abs. value of difference between mean endogenous PCo1-3 and LS mean PCo1-3 for each age group 
###    at each site, fit linear functions of age, save Betas (do for both leaves and roots)
#### Test prediction: Roots get less similar to roots in GRH potting soil with age
### 8. Calculate abs. value of difference between mean potting soil PCo1 and LS mean PCo1 for each age group
###    at each site, fit linear functions of age, save Betas (no data for leaves, so roots only)
### 9. Repeat 5-8 for 1000 datasets with Age permuted; compare observed Betas to distributions

####### Define zero-tolerant function for calculating geometric means (credit: P.J. McMurdie) #######
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
####### variance stabilizing transformation of FULL dataset including greenhouse samples #######
## Apply variance stabilizing transformation
fullEco.vst<-fullEco # copy Phyloseq object
fullEco.dds<-phyloseq_to_deseq2(fullEco.vst,~Type+Site) # Change to DESeq2 object-- can't include Type*Site because there are no leaves from the greenhouse
fullEco.dds = estimateSizeFactors(fullEco.dds, geoMeans = apply(counts(fullEco.dds), 1, gm_mean)) # Estimate size factors using geometric means 
fullEco.rld <- DESeq2::varianceStabilizingTransformation(fullEco.dds, blind = FALSE) # Apply variance stabilizing transformation
fullEco.vstMat <- GenomicRanges::assay(fullEco.rld) # Extract transformed OTU table
otu_table(fullEco.vst) <- otu_table(fullEco.vstMat, taxa_are_rows = TRUE) # Put transformed OTU table back into Phyloseq object
rm(fullEco.dds,fullEco.rld,fullEco.vstMat) # Clean up
save(fullEco.vst,file='intermediate_data/Phyloseq_fullEco_vst.RData')

####### Ordinate on weighted UniFrac distance ####### 

registerDoParallel(cores=4)

fullEco.vst.pos<-transform_sample_counts(fullEco.vst,function(x) x<-ifelse(x<0,0,x)) # round up negative values for compatibility with ordination methods
rootsEco.vst.pos<-subset_samples(fullEco.vst.pos,Type!='leaf') # only ordinate roots and soil

wUF.rootsEco.vst.pos<-UniFrac(rootsEco.vst.pos,weighted=TRUE,parallel=TRUE) # calculate weighted UniFrac distances
save(wUF.rootsEco.vst.pos,file='succession/wUF_allrootsoilEco_vst_pos.RData')
PCoA.wUF.rootsEco.vst<-capscale(wUF.rootsEco.vst.pos~1) # ordinate based on weighted UniFrac distances

## View  proportion variance explained by each PCo 
summary(PCoA.wUF.rootsEco.vst) %>% head()

## Plot ordination of root and soil samples from field and greenhouse
pdf(file="plots/Fig_S10a_Ordination_wUF_rootandsoil.pdf",width=9,height=9)
subset_samples(rootsEco.vst.pos,Treatment%in%c('field','GRHsoil') & !(Site%in%c('Mil','Par'))) %>%
plot_ordination(.,PCoA.wUF.rootsEco.vst,color='Site',shape='Type')+
  geom_point(size=3)+
  theme_classic()+
  xlab("PCo1 [27.4%]")+ylab("PCo2 [13.9%]")+
  theme(axis.title.x=element_text(size=26,face="bold"),axis.text.x=element_text(size=20))+
  theme(axis.title.y=element_text(size=26,face="bold"),axis.text.y=element_text(size=20))+
  theme(legend.title=element_text(size=26,face='bold'),legend.text=element_text(size=22,face='bold'))+
  #scale_color_manual(values=c('forest green',sitePalette[1:2],'grey','purple',sitePalette[3]))+
  scale_color_manual(values=c('forest green',sitePalette),labels=c('Greenhouse','Jam','Mah','Sil'))+
  scale_shape_manual(values=c(16:17),name='')
dev.off()

## Save PCoA axes for use in permutation tests 
roots.smd<-as(sample_data(rootsEco.vst.pos),'data.frame')
roots.smd<-scores(PCoA.wUF.rootsEco.vst,display="sites",choices=c(1:3)) %>% as.data.frame() %>% 
  rename(wUF1=MDS1,wUF2=MDS2,wUF3=MDS3) %>% 
  merge(roots.smd,.,by.x='SampleID',by.y='row.names')

####### Separate out roots growing in normal potting soil #######
grhsoil.smd<-filter(roots.smd,Type=='root',Treatment=='GRHsoil')
roots.smd<-filter(roots.smd,Treatment!='GRHsoil') # remove from field plants

####### Separate out roots of endogenous plants #######
endog.smd<-filter(roots.smd,Type=='root',Age=='endog')
roots.smd<-filter(roots.smd,Age!='endog') # remove from experimental plants

# Get mean PCoA coordinates for endogenous roots at each Site
endog.summ<-summarize(group_by(endog.smd,Site),meanwUF1=mean(wUF1),meanwUF2=mean(wUF2),meanwUF3=mean(wUF3))

####### Only keep roots (exclude soils) and exclude field sites with low survival #######
roots.smd<-filter(roots.smd,Type=='root',!(Site%in%c('Mil','Par')))
# these are now the same samples as in 'root3.smd' in main code, i.e.,
## all roots from the 3 main field sites with good survival

####### Permutations: weighted UniFrac: leave all variables the same except permute Age #######
## IF succession is still happening, we expect: 
### 1. Sites get more distinct with age (in roots) (i.e., variance among sites increases with age)
### 2. Roots but not leaves get more similar to endogenous plants with age
### 3. Roots get less similar to GRH potting soil

set.seed(777)
nperm<-1000 # set number of permutations: higher than 1000 because some might fail to converge-- will then choose the first 1000 successful permutations from the set.
perm.roots.smd<-roots.smd # copy data frame to be permuted
i<-1
while(i<nperm) {
if (i==1){betas.wUF<-data.frame("Simulation"=integer(),
  "betaVarSites.wUF1"=numeric(),"betaVarSites.wUF2"=numeric(),"betaVarSites.wUF3"=numeric(),
  "betaDistGRH.wUF1"=numeric(),"betaDistGRH.wUF2"=numeric(),"betaDistGRH.wUF3"=numeric(),
  "betaDistEndog.wUF1"=numeric(),"betaDistEndog.wUF2"=numeric(),"betaDistEndog.wUF3"=numeric())
              } # initialize/clear output data frame for first iteration
  print(paste0("Beginning simulation number ",as.character(i)))
  print(paste0("Fitting LMMs... iteration ",as.character(i)))
  perm.roots.smd$Age<-sample(perm.roots.smd$Age,length(perm.roots.smd$Age),replace=FALSE) # permute Ages only
  perm.root3.wUF1.SxA<-try(lmer(wUF1~Genotype*Site + Site*Age + Site*Harvested + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=perm.roots.smd,REML=TRUE)) # fit wUF1 model
  perm.root3.wUF2.SxA<-try(lmer(wUF2~Genotype*Site + Site*Age + Site*Harvested + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=perm.roots.smd,REML=TRUE)) # fit wUF2 model
  perm.root3.wUF3.SxA<-try(lmer(wUF3~Genotype*Site + Site*Age + Site*Harvested + (1|Site:Block) + (1|Genotype:Line)+(1|newPlate)+logObs,data=perm.roots.smd,REML=TRUE)) # fit wUF3 model
  
  print(paste0("Extracting LS means..... iteration ",as.character(i)))
  # Try to extract LS means. If model didn't converge, the LSM object will be of class 'try-error'
  LSM.perm.wUF1<-try(as.data.frame(lmerTest::lsmeans(perm.root3.wUF1.SxA)[1])) 
  LSM.perm.wUF2<-try(as.data.frame(lmerTest::lsmeans(perm.root3.wUF2.SxA)[1]))
  LSM.perm.wUF3<-try(as.data.frame(lmerTest::lsmeans(perm.root3.wUF3.SxA)[1]))

  ## Determine whether all 3 models converged for this permutation:
  status<-(ifelse(class(LSM.perm.wUF1)=='try-error' | class(LSM.perm.wUF2)=='try-error' | class(LSM.perm.wUF3)=='try-error', 
             'fail', 'success'))
  print(paste0("Permutation ",i," status = ",status))
  
  ### Only proceed if all 3 models converged: ###
  if (status=='success') { 
    # Combine all LS means into one data frame
    LSM.perm.all<-rbind(mutate(LSM.perm.wUF1,Response="wUF1",Term=rownames(LSM.perm.wUF1)),
                       mutate(LSM.perm.wUF2,Response="wUF2",Term=rownames(LSM.perm.wUF2)),
                       mutate(LSM.perm.wUF3,Response="wUF3",Term=rownames(LSM.perm.wUF3)))
    colnames(LSM.perm.all)<-gsub('lsmeans.table.','',colnames(LSM.perm.all))

    # For testing prediction that variance among sites increases with age (sites getting more distinct from each other as common greenhouse-based component decreases)
    print(paste0("Calculating variance among sites....... iteration ",as.character(i)))
    varwUF.perm<-filter(LSM.perm.all,grepl("Site:Age",Term)) %>% group_by(.,Age,Response) %>% 
      summarize(.,VarianceAmongSites=var(Estimate))
  
    betaVarSites.wUF1<-coef(lm(VarianceAmongSites~as.numeric(Age),data=filter(varwUF.perm,Response=='wUF1')))[2]
    betaVarSites.wUF2<-coef(lm(VarianceAmongSites~as.numeric(Age),data=filter(varwUF.perm,Response=='wUF2')))[2]
    betaVarSites.wUF3<-coef(lm(VarianceAmongSites~as.numeric(Age),data=filter(varwUF.perm,Response=='wUF3')))[2]
    
    # For testing prediction that roots in the field become less similar to roots in potting soil with age
    print(paste0("Calculating distance from roots in potting soil.... iteration ",as.character(i)))
    distGRH.perm<-filter(LSM.perm.all,grepl("Site:Age",Term)) %>%
      mutate(distGRH=ifelse(Response=='wUF1',abs(mean(grhsoil.smd$wUF1)-Estimate), # calculate absolute value of distance from mean potting soil root
                            ifelse(Response=='wUF2',abs(mean(grhsoil.smd$wUF2)-Estimate), # do separately for each site to allow sites to move in different directions
                                   abs(mean(grhsoil.smd$wUF3)-Estimate))))
    betaDistGRH.wUF1<-coef(lm(distGRH~as.numeric(Age),data=filter(distGRH.perm,Response=='wUF1')))[2]
    betaDistGRH.wUF2<-coef(lm(distGRH~as.numeric(Age),data=filter(distGRH.perm,Response=='wUF2')))[2]
    betaDistGRH.wUF3<-coef(lm(distGRH~as.numeric(Age),data=filter(distGRH.perm,Response=='wUF3')))[2]
    
    # For testing prediction that roots in each site get more similar to local endogenous roots with age 
    print(paste0("Calculating distance from endogenous plants at each site.... iteration ",as.character(i)))
    distEndog.perm<-filter(LSM.perm.all,grepl("Site:Age",Term)) %>% merge(endog.summ,by='Site') %>%
      mutate(distEndog=ifelse(Response=='wUF1',abs(meanwUF1-Estimate), # calculate absolute value of distance from endogenous mean at each site
                              ifelse(Response=='wUF2',abs(meanwUF2-Estimate), # do separately for each site to allow sites to move in different directions
                                     abs(meanwUF3-Estimate))))
    betaDistEndog.wUF1<-coef(lm(distEndog~as.numeric(Age),data=filter(distEndog.perm,Response=='wUF1')))[2]
    betaDistEndog.wUF2<-coef(lm(distEndog~as.numeric(Age),data=filter(distEndog.perm,Response=='wUF2')))[2]
    betaDistEndog.wUF3<-coef(lm(distEndog~as.numeric(Age),data=filter(distEndog.perm,Response=='wUF3')))[2]
    
    ### Add to dataframe of betas:
    betas.wUF<-rbind(betas.wUF,
                 data.frame("Simulation"=i,
                            "betaVarSites.wUF1"=betaVarSites.wUF1,"betaVarSites.wUF2"=betaVarSites.wUF2,"betaVarSites.wUF3"=betaVarSites.wUF3,
                            "betaDistGRH.wUF1"=betaDistGRH.wUF1,"betaDistGRH.wUF2"=betaDistGRH.wUF2,"betaDistGRH.wUF3"=betaDistGRH.wUF3,
                            "betaDistEndog.wUF1"=betaDistEndog.wUF1,"betaDistEndog.wUF2"=betaDistEndog.wUF2,"betaDistEndog.wUF3"=betaDistEndog.wUF3)) # store simulation results in growing dataframe
  rm(LSM.perm.all) # clean up LS means
  rm(varwUF.perm,betaVarSites.wUF1,betaVarSites.wUF2,betaVarSites.wUF3) # clean up variance estimates 
  rm(distGRH.perm,betaDistGRH.wUF1,betaDistGRH.wUF2,betaDistGRH.wUF3) # clean up estimates of distance from greenhouse roots
  rm(distEndog.perm,betaDistEndog.wUF1,betaDistEndog.wUF2,betaDistEndog.wUF3) # clean up estimates of distance from local endogenous roots
  i <- i+1 # advance the counter only if the simulation was successful (i.e., the models converged)
  }
  rm(perm.root3.wUF1.SxA,perm.root3.wUF2.SxA,perm.root3.wUF3.SxA) # clean up models
  rm(LSM.perm.wUF1,LSM.perm.wUF2,LSM.perm.wUF3) # clean up LSM objects
}

save(betas.wUF,file="succession/betas_wUF.RData")

####### Weighted UniFrac: Fit standard model of PCo1-3 from 3 sites and extract Age*Site LS means #######
obs1<-lmer(wUF1~Genotype*Site+Age*Site+Harvested*Site+(1|Site:Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=roots.smd,REML=TRUE)
qqnorm(resid(obs1));qqline(resid(obs1))
obs2<-lmer(wUF2~Genotype*Site+Age*Site+Harvested*Site+(1|Site:Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=roots.smd,REML=TRUE)
qqnorm(resid(obs2));qqline(resid(obs2))
obs3<-lmer(wUF3~Genotype*Site+Age*Site+Harvested*Site+(1|Site:Block)+(1|Genotype:Line)+(1|newPlate)+logObs,data=roots.smd,REML=TRUE)
qqnorm(resid(obs3));qqline(resid(obs3))

anova(obs1)
anova(obs2)
anova(obs3)
# All 3 axes show significant Site:Age interaction

LSM.obs.wUF1<-as.data.frame(lmerTest::lsmeans(obs1)[1])
LSM.obs.wUF2<-as.data.frame(lmerTest::lsmeans(obs2)[1])
LSM.obs.wUF3<-as.data.frame(lmerTest::lsmeans(obs3)[1])
LSM.obs.all<-rbind(mutate(LSM.obs.wUF1,Response="wUF1",Term=rownames(LSM.obs.wUF1)),
                   mutate(LSM.obs.wUF2,Response="wUF2",Term=rownames(LSM.obs.wUF2)),
                   mutate(LSM.obs.wUF3,Response="wUF3",Term=rownames(LSM.obs.wUF3)))
colnames(LSM.obs.all)<-gsub('lsmeans.table.','',colnames(LSM.obs.all)) # shorten column names
rm(LSM.obs.wUF1,LSM.obs.wUF2,LSM.obs.wUF3)

####### Weighted UniFrac: Calculate observed values: linear change in variance among sites with age #######
# Calculate variance among sites for each PCoA axis for each age group, using LS means
varwUF<-filter(LSM.obs.all,grepl("Site:Age",Term)) %>% group_by(.,Age,Response) %>% 
  summarize(.,VarianceAmongSites=var(Estimate))

### Positive coefficients are consistent with ongoing succession: 
trueBetaVarSites.wUF1<-coef(lm(VarianceAmongSites~as.numeric(Age),data=filter(varwUF,Response=='wUF1')))[2] # positive coefficient
trueBetaVarSites.wUF2<-coef(lm(VarianceAmongSites~as.numeric(Age),data=filter(varwUF,Response=='wUF2')))[2] # negative coefficient
trueBetaVarSites.wUF3<-coef(lm(VarianceAmongSites~as.numeric(Age),data=filter(varwUF,Response=='wUF3')))[2] # negative coefficient

####### Weighted UniFrac: Calculate observed values: linear change in distance from potting soil roots with age #######
distGRH<-filter(LSM.obs.all,grepl("Site:Age",Term)) %>%
  mutate(distGRH=ifelse(Response=='wUF1',abs(mean(grhsoil.smd$wUF1)-Estimate), # calculate absolute value of distance from mean potting soil root
                        ifelse(Response=='wUF2',abs(mean(grhsoil.smd$wUF2)-Estimate), # do separately for each site to allow sites to move in different directions
                               abs(mean(grhsoil.smd$wUF3)-Estimate))))

### Positive coefficients are consistent with ongoing succession: 
trueBetaDistGRH.wUF1<-coef(lm(distGRH~as.numeric(Age),data=filter(distGRH,Response=='wUF1')))[2] # positive
trueBetaDistGRH.wUF2<-coef(lm(distGRH~as.numeric(Age),data=filter(distGRH,Response=='wUF2')))[2] # negative
trueBetaDistGRH.wUF3<-coef(lm(distGRH~as.numeric(Age),data=filter(distGRH,Response=='wUF3')))[2] # positive

####### Weighted UniFrac: Calculate observed values: linear change in distance to endogenous plants at each site, with age #######
distEndog<-filter(LSM.obs.all,grepl("Site:Age",Term)) %>% merge(endog.summ,by='Site') %>%
  mutate(distEndog=ifelse(Response=='wUF1',abs(meanwUF1-Estimate), # calculate absolute value of distance from endogenous mean at each site
                          ifelse(Response=='wUF2',abs(meanwUF2-Estimate), # do separately for each site to allow sites to move in different directions
                                 abs(meanwUF3-Estimate))))

### Negative coefficients are consistent with ongoing succession: 
trueBetaDistEndog.wUF1<-coef(lm(distEndog~as.numeric(Age),data=filter(distEndog,Response=='wUF1')))[2] # negative
trueBetaDistEndog.wUF2<-coef(lm(distEndog~as.numeric(Age),data=filter(distEndog,Response=='wUF2')))[2] # negative
trueBetaDistEndog.wUF3<-coef(lm(distEndog~as.numeric(Age),data=filter(distEndog,Response=='wUF3')))[2] # positive

####### Histograms: Test predictions: increasing variance among sites  #######
## Prediction #1: if succession is ongoing, we expect: 
##### 1. The average root communities at the 3 sites get more distinct with age (because they started in a common potting-soil inoculum)
pdf(file="plots/Fig_S11a_needsPvals.pdf",width=10,height=9)
select(betas.wUF,starts_with('betaVarSites')) %>% 
  melt(measure.vars=colnames(.)[1:3], variable.name="Axis",value.name="Position") %>% 
  mutate(Axis=gsub("betaVarSites.wUF","PCo",Axis)) %>% 
  mutate(TrueBeta=ifelse(Axis=='PCo1',trueBetaVarSites.wUF1,
                         ifelse(Axis=='PCo2',trueBetaVarSites.wUF2,trueBetaVarSites.wUF3))) %>% 
  ggplot(.,aes(x=Position))+
  facet_grid(Axis~.,scales="free")+
  geom_histogram(binwidth=0.001,fill='grey',color='black')+
  geom_vline(aes(xintercept=TrueBeta),color='red',size=1.5,lty=2)+
  xlab("(i) Slope: Variance among sites ~ age")+
  theme_classic()+
  theme(strip.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))+
  theme(axis.title.x=element_text(size=30,face="bold"),axis.text.x=element_text(size=20))+
  theme(axis.title.y=element_text(size=26,face="bold"),axis.text.y=element_text(size=20))
dev.off()
####### Histograms: Test predictions: decreasing dissimilarity to local endogenous plants #######
## Prediction #2: 
## IF succession is still happening, we expect: roots get more similar to endogenous plants with age
pdf(file="plots/Fig_S11b_needsPvals.pdf",width=10,height=9)
select(betas.wUF,starts_with('betaDistEndog')) %>% 
  melt(measure.vars=colnames(.)[1:3], variable.name="Axis",value.name="Position") %>% 
  mutate(Axis=gsub("betaDistEndog.wUF","PCo",Axis)) %>% 
  mutate(TrueBeta=ifelse(Axis=='PCo1',trueBetaDistEndog.wUF1,
                         ifelse(Axis=='PCo2',trueBetaDistEndog.wUF2,trueBetaDistEndog.wUF3))) %>% 
  ggplot(.,aes(x=Position))+
  facet_grid(Axis~.,scales="free")+
  geom_histogram(binwidth=0.005,fill='grey',color='black')+
  geom_vline(aes(xintercept=TrueBeta),color='red',size=1.5,lty=2)+
  xlab("(ii) Slope: Distance from endogenous roots ~ age")+
  theme_classic()+
  theme(strip.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))+
  theme(axis.title.x=element_text(size=30,face="bold"),axis.text.x=element_text(size=20))+
  theme(axis.title.y=element_text(size=26,face="bold"),axis.text.y=element_text(size=20))
dev.off()
####### Histograms: Test predictions: decreasing dissimilarity to roots in potting-soil #######
## Prediction #3: 
## IF succession is still happening, we expect: roots get less similar to roots in potting-soil as plants age
pdf(file="plots/Fig_S11c_needsPvals.pdf",width=10,height=9)
select(betas.wUF,starts_with('betaDistGRH')) %>% 
  melt(measure.vars=colnames(.)[1:3], variable.name="Axis",value.name="Position") %>% 
  mutate(Axis=gsub("betaDistGRH.wUF","PCo",Axis)) %>% 
  mutate(TrueBeta=ifelse(Axis=='PCo1',trueBetaDistGRH.wUF1,
                     ifelse(Axis=='PCo2',trueBetaDistGRH.wUF2,trueBetaDistGRH.wUF3))) %>% 
ggplot(.,aes(x=Position))+
  facet_grid(Axis~.,scales="free")+
  geom_histogram(binwidth=0.005,fill='grey',color='black')+
  geom_vline(aes(xintercept=TrueBeta),color='red',size=1.5,lty=2)+
  xlab("(iii) Slope: Distance from greenhouse roots ~ age")+
  theme_classic()+
  theme(strip.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))+
  theme(axis.title.x=element_text(size=30,face="bold"),axis.text.x=element_text(size=20))+
  theme(axis.title.y=element_text(size=26,face="bold"),axis.text.y=element_text(size=20))
dev.off()

####### P-values: Testing predictions #######
nsim<-max(betas.wUF$Simulation)
#### First test the expectation that the observed (true) change in variance among sites with age is 
#### equally or more positive than in 999 simulated datasets with permuted age.
#### P-value = # of simulated slopes > true slope, / nsim (= number of permutations)
filter(betas.wUF,betaVarSites.wUF1>=trueBetaVarSites.wUF1) %>% select(betaVarSites.wUF1) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo1  P=0
filter(betas.wUF,betaVarSites.wUF2>=trueBetaVarSites.wUF2) %>% select(betaVarSites.wUF2) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo2  P=1
filter(betas.wUF,betaVarSites.wUF3>=trueBetaVarSites.wUF3) %>% select(betaVarSites.wUF3) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo3  P=0.71

#### Second, test the expectation that the observed (true) change in dissimilarity to endogenous roots with age is 
#### equally or more NEGATIVE than in 999 simulated datasets with permuted age (i.e., if succession is happening,
#### age will bring the root communities more similar to those in endogenous plants)
#### P-value = # of simulated slopes < true slope, / 999 (= number of permutations)
filter(betas.wUF,betaDistEndog.wUF1<trueBetaDistEndog.wUF1) %>% select(betaDistEndog.wUF1) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo1  P=0
filter(betas.wUF,betaDistEndog.wUF2<trueBetaDistEndog.wUF2) %>% select(betaDistEndog.wUF2) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo2  P=0
filter(betas.wUF,betaDistEndog.wUF3<trueBetaDistEndog.wUF3) %>% select(betaDistEndog.wUF3) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo3  P=0.84

#### Third, test the expectation that the observed (true) change in dissimilarity to greenhouse roots with age is 
#### equally or more POSITIVE than in 999 simulated datasets with permuted age (i.e., if succession is happening,
#### age will bring the root communities LESS similar to those in greenhouse/potting-soil plants)
#### P-value = # of simulated slopes > true slope, / 999 (= number of permutations)
filter(betas.wUF,betaDistGRH.wUF1>=trueBetaDistGRH.wUF1) %>% select(betaDistGRH.wUF1) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo1  P=0
filter(betas.wUF,betaDistGRH.wUF2>=trueBetaDistGRH.wUF2) %>% select(betaDistGRH.wUF2) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo2  P=1
filter(betas.wUF,betaDistGRH.wUF3>=trueBetaDistGRH.wUF3) %>% select(betaDistGRH.wUF3) %>%
  dim() %>% as.data.frame() %>% slice(1)/nsim # PCo3  P=0.003

####### Save image #######
save.image(paste0("succession/succession_image",date(),".RData"))
