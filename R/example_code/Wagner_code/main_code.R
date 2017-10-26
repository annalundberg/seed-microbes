### Author: Maggie R. Wagner
### maggie.r.wagner@gmail.com
####### Move all data files into raw_data subdirectory #######
system('mkdir raw_data')
system('mv contaminants.fasta raw_data')
system('mv Eco_Field_copynumest_forR.txt raw_data')
system('mv Ecotypes_field_glucosinolates.txt raw_data')
system('mv OTUrepSeqs97.fa raw_data')
system('mv otuTable97.txt raw_data')
system('mv otuTable99.txt raw_data')
system('mv phylogeny.tre raw_data')
system('mv plant_key.txt raw_data')
system('mv site_coords_Ecotypes.txt raw_data')
system('mv SMD.txt raw_data')
system('mv soildata.txt raw_data')
system('mv taxAssignments97.txt raw_data')
system('mv taxAssignments99.txt raw_data')

####### Separate code sections into subdirectories #######
system('mkdir foldchange')
system('mkdir GLS')
system('mkdir heritability')
system('mkdir higher_tax_levels')
system('mkdir LMMs')
system('mkdir otu99pct')
system('mkdir rel_abund')
system('mkdir succession')
system('mkdir uUF')
# Move scripts into respective folders:
system('mv foldchange.R foldchange')
system('mv GLS.R GLS')
system('mv heritability.R heritability')
system('mv countmodels_vstLMM.R heritability')
system('mv populate_vstLMM.R heritability')
system('mv higher_tax_levels.R higher_tax_levels')
system('mv LMMs.R LMMs')
system('mv otu99pct.R otu99pct')
system('mv rel_abund.R rel_abund')
system('mv succession.R succession')
system('mv uUF.R uUF')
####### Make subdirectories for plots, tables, and intermediate data files #######
system('mkdir plots')
system('mkdir tables')
system('mkdir intermediate_data')

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
[49] acepack_1.3-3.3     "
####### Load source file ######
source('ecotypes_source.R')
####### Load individual components of Phyloseq object #######
# Import raw OTU table (not copy number corrected yet)
fullOTU<-as.data.frame(read.table('raw_data/otuTable97.txt',sep='\t',header=TRUE))
rownames(fullOTU)<-paste('OTU_',as.character(fullOTU[,1]),sep='') # store OTU_IDs as row names
fullOTU<-fullOTU[,-1] # remove 'OTU_ID' column
fullOTU<-data.matrix(fullOTU) # convert to data matrix format

## load sample metadata
fullSMD<-as.data.frame(read.table('raw_data/SMD.txt',sep='\t',header=TRUE))
# drop columns that won't be used:
fullSMD<-fullSMD[,!names(fullSMD) %in% c("oldPlate")]
# clean up data: (change Year, Age to factors and place endogenous plants into Ecotypes experiment)
endog<-filter(fullSMD,Genotype=='endog')
endog$Age<-as.factor('endog')
endog$Cohort<-as.factor('endog')

fullSMD<-filter(fullSMD,Genotype!='endog')%>%rbind(.,endog)
fullSMD$Age<-factor(fullSMD$Age)
fullSMD$Harvested<-factor(fullSMD$Harvested)
fullSMD$Experiment<-plyr::mapvalues(fullSMD$Experiment,c('endog'),c('ecotypes'))

rownames(fullSMD)<-fullSMD$SampleID
## load taxonomy assignments
taxfile<-as.matrix(read.table('raw_data/taxAssignments97.txt',sep='\t',header=TRUE))
dimnames(taxfile)[[1]]<-taxfile[,1]  # make OTU_IDs the row names
taxfile<-taxfile[,-c(1,2,8)]  # get rid of taxonomy, OTU.ID, Confidence columns
taxfile<- as.data.frame(taxfile)
# Fix problem noticed during downstream data analysis:
filter(taxfile,Family=='Chromatiaceae') # Family Chromatiaceae classified into 2 different orders: Chromatiales (correct) and Alteromonadales (incorrect)
badChromatiaceae<-subset(taxfile,Family=='Chromatiaceae' & Order=='Alteromonadales')
# Best BLAST matches for 'Chromatiaceae-Alteromonadales':
# OTU_4719 = Rheinheimera sp. = Chromatiaceae, Chromatiales
# OTU_29684 = Alishewanella sp. = Alteromonadaceae, Alteromonadales
# OTU_22285 = Rheinheimera sp. = Chromatiaceae, Chromatiales

# Fix manually: 
badChromatiaceae$Order<-ifelse(rownames(badChromatiaceae)%in%c('OTU_4719','OTU_22285'),'Chromatiales','Alteromonadales')
badChromatiaceae$Family<-ifelse(rownames(badChromatiaceae)=='OTU_29684','Alteromonadaceae','Chromatiaceae')
taxfile<-subset(taxfile,!(Family=='Chromatiaceae' & Order=='Alteromonadales')) %>%
  rbind(.,badChromatiaceae) 
taxfile<-as.matrix(taxfile)
rm(badChromatiaceae)

## load phylogeny
phyfile<-read.tree(file="raw_data/phylogeny.tre")

####### Load reference sequences and contaminants ######
## load reference sequences
refseq<-readDNAStringSet("raw_data//OTUrepSeqs97.fa",format="fasta")
contams<-readDNAStringSet("raw_data//contaminants.fasta",format="fasta")
# Extract the OTU_IDs of the contaminants:
contams<-names(contams) %>% strsplit(.," ") %>% unlist() # extract sequence names
contams<-contams[grep("OTU_",contams)] # only keep the OTU_IDs

####### Generate master phyloseq object with ALL samples & data collected in 2011-2012 #######
full97<-phyloseq(otu_table(fullOTU,taxa_are_rows=TRUE),sample_data(fullSMD),tax_table(taxfile),phyfile,refseq)
rm(fullOTU,taxfile,fullSMD,phyfile) # clean up individual components
## Get sequencing depth for each sample (before copy number correction or removing bad OTUs)
sample_data(full97)$HQReads<-sample_sums(full97) # add up high-quality reads in each sample
sample_data(full97)$logHQReads<-log(sample_data(full97)$HQReads) # store natural log of read counts
####### Save only the samples involved in the "Ecotypes" experiment #######
fullEco<-subset_samples(full97,Analysis=='Ecotypes')
save(fullEco,file="intermediate_data/fullEco_noCNCnoQC.RData")
####### Remove host sequences/bad OTUs and contaminants #######
fullEco.nobadOTUs<-subset_taxa(fullEco,Kingdom!='Unassigned' & Class!='Chloroplast' & Family !='mitochondria')
fullEco.nobadOTUs<-prune_taxa(setdiff(taxa_names(fullEco.nobadOTUs),contams),fullEco.nobadOTUs) # remove contaminant OTUs
  
sample_data(fullEco.nobadOTUs)$UsableReads<-sample_sums(fullEco.nobadOTUs)
####### Remove samples with <800 usable reads #######
fullEco.nobadOTUs.highcoverage<-subset_samples(fullEco.nobadOTUs,UsableReads>=800)
save(fullEco.nobadOTUs.highcoverage,file="intermediate_data/fullEco_preThreshold_preCNC.RData")

####### Subset out samples from the Field experiment only (separate out greenhouse samples) #######
fieldEco.nobadOTUs.highcoverage<-subset_samples(fullEco.nobadOTUs.highcoverage,Site!='Duke')
save(fieldEco.nobadOTUs.highcoverage,file="intermediate_data/fieldEco_preThreshold_preCNC.RData")

####### Standard thresholding: 5x25 (throw out "non-reproducible" OTUs) #######
threshold1<-kOverA(5,A=25) # set threshold values (require k samples with A reads)
fullEco.nobadOTUs.highcoverage.thresholded<-filter_taxa(fullEco.nobadOTUs.highcoverage,threshold1,TRUE)
fieldEco.nobadOTUs.highcoverage.thresholded<-filter_taxa(fieldEco.nobadOTUs.highcoverage,threshold1,TRUE)
####### Copy number correction #######
## Load copy number estimates:
CNCest<-as.data.frame(read.table('raw_data/Eco_Field_copynumest_forR.txt',sep='\t',header=TRUE))

# do for fullEco dataset:
fullOTU.CNC<-as.data.frame(as(otu_table(fullEco.nobadOTUs.highcoverage.thresholded),'matrix')) # copy OTU table from thresholded Phyloseq object
fullOTU.CNC<-merge(CNCest,fullOTU.CNC,by.y='row.names',by.x='OTU_ID') # merge copy number estimates with OTU table
rownames(fullOTU.CNC)<-fullOTU.CNC$OTU_ID # restore OTU_IDs as row names
fullOTU.CNC<-fullOTU.CNC[,3:dim(fullOTU.CNC)[2]]/fullOTU.CNC[,2] # divide OTU counts by corresponding copy number estimate (Col 1 = OTU_ID, Col 2 = copy number estimates)
fullOTU.CNC<-data.matrix(ceiling(fullOTU.CNC)) # round OTU counts up to nearest integer

# do for fieldEco dataset:
fieldOTU.CNC<-as.data.frame(as(otu_table(fieldEco.nobadOTUs.highcoverage.thresholded),'matrix')) # copy OTU table from thresholded Phyloseq object
fieldOTU.CNC<-merge(CNCest,fieldOTU.CNC,by.y='row.names',by.x='OTU_ID') # merge copy number estimates with OTU table
rownames(fieldOTU.CNC)<-fieldOTU.CNC$OTU_ID # restore OTU_IDs as row names
fieldOTU.CNC<-fieldOTU.CNC[,3:dim(fieldOTU.CNC)[2]]/fieldOTU.CNC[,2] # divide OTU counts by corresponding copy number estimate (Col 1 = OTU_ID, Col 2 = copy number estimates)
fieldOTU.CNC<-data.matrix(ceiling(fieldOTU.CNC)) # round OTU counts up to nearest integer

## Replace OTU table of master Phyloseq objects with copy number corrected versions:
# do for fullEco dataset:
fullEco.nobadOTUs.highcoverage.thresholded.CNC<-fullEco.nobadOTUs.highcoverage.thresholded # copy Phyloseq object
otu_table(fullEco.nobadOTUs.highcoverage.thresholded.CNC)<-otu_table(fullOTU.CNC,taxa_are_rows=TRUE) # replace OTU table with copy number-corrected version

# do for fieldEco dataset:
fieldEco.nobadOTUs.highcoverage.thresholded.CNC<-fieldEco.nobadOTUs.highcoverage.thresholded # copy Phyloseq object
otu_table(fieldEco.nobadOTUs.highcoverage.thresholded.CNC)<-otu_table(fieldOTU.CNC,taxa_are_rows=TRUE) # replace OTU table with copy number-corrected version

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(fieldEco.nobadOTUs.highcoverage.thresholded.CNC)$Obs<-sample_sums(fieldEco.nobadOTUs.highcoverage.thresholded.CNC) # add up "good" OTU observations in each sample
sample_data(fieldEco.nobadOTUs.highcoverage.thresholded.CNC)$logObs<-log(sample_data(fieldEco.nobadOTUs.highcoverage.thresholded.CNC)$Obs) # store natural log of "good" OTU observations

sample_data(fullEco.nobadOTUs.highcoverage.thresholded.CNC)$Obs<-sample_sums(fullEco.nobadOTUs.highcoverage.thresholded.CNC) # add up "good" OTU observations in each sample
sample_data(fullEco.nobadOTUs.highcoverage.thresholded.CNC)$logObs<-log(sample_data(fullEco.nobadOTUs.highcoverage.thresholded.CNC)$Obs) # store natural log of "good" OTU observations

# how many observations?
mean(sample_data(fieldEco.nobadOTUs.highcoverage.thresholded.CNC)$Obs)
sd(sample_data(fieldEco.nobadOTUs.highcoverage.thresholded.CNC)$Obs)

# how many OTUs?
ntaxa(fieldEco.nobadOTUs.highcoverage.thresholded.CNC)

####### Save fully cleaned datasets and clean up Phyloseq objects in various other stages of QC ########
save(fieldEco.nobadOTUs.highcoverage.thresholded.CNC,file="intermediate_data/fieldEco_cleaned_CNC.RData")
save(fullEco.nobadOTUs.highcoverage.thresholded.CNC,file="intermediate_data/fullEco_cleaned_CNC.RData")
rm(phyfile,full97,taxfile,fullOTU,fullSMD,threshold1)
rm(CNCest,fieldEco.nobadOTUs.highcoverage,fieldEco.nobadOTUs.highcoverage.thresholded,fieldOTU.CNC) # clean up- all data are now stored in fieldEco...CNC Phyloseq object
rm(fullEco,fullEco.nobadOTUs,fullEco.nobadOTUs.highcoverage,fullEco.nobadOTUs.highcoverage.thresholded,fullOTU.CNC) # clean up- all data are now stored in fullEco...CNC Phyloseq object

####### Consolidate data at higher taxonomic levels #######
# execute './higher_tax_levels/higher_tax_levels.R'

####### Calculate relative abundances of OTUs, families, ..., phyla #######
# execute './rel_abund/rel_abund.R'

####### Analyses to disentangle post-transplant succession from plant age #######
source('./succession/succession.R')

####### Remove "non-endogenous" OTUs from the field dataset (see succession.R) #######
load("succession/non_endogenous_OTUs.RData")
ntaxa(fieldEco.nobadOTUs.highcoverage.thresholded.CNC) #3718 OTUs to start
fieldEco.nobadOTUs.highcoverage.thresholded.CNC<-prune_taxa(
  setdiff(taxa_names(fieldEco.nobadOTUs.highcoverage.thresholded.CNC),nonEndog.OTUs),
  fieldEco.nobadOTUs.highcoverage.thresholded.CNC)
ntaxa(fieldEco.nobadOTUs.highcoverage.thresholded.CNC) #3689 remain (29 were removed)

####### Store OTU taxonomic table for easy merging in future #######
OTUtax<-as(tax_table(fieldEco.nobadOTUs.highcoverage.thresholded.CNC),'matrix') %>%
  as.data.frame %>% mutate(OTU_ID=row.names(.))
save(OTUtax,file='intermediate_data/OTU_tax_table.RData')

####### # variance stabilizing transformation of entire field dataset: #######
## Apply variance stabilizing transformation
## Warning: this transformation takes multiple days to complete. Recommended to execute on a HPC cluster
fieldEco.vst<-fieldEco.nobadOTUs.highcoverage.thresholded.CNC # copy Phyloseq object
fieldEco.dds<-phyloseq_to_deseq2(fieldEco.vst,~Type+Site+Type*Site) # Change to DESeq2 object
fieldEco.dds = estimateSizeFactors(fieldEco.dds, geoMeans = apply(counts(fieldEco.dds), 1, gm_mean)) # Estimate size factors using geometric means 
fieldEco.rld <- DESeq2::varianceStabilizingTransformation(fieldEco.dds, blind = FALSE) # Apply variance stabilizing transformation
fieldEco.vstMat <- GenomicRanges::assay(fieldEco.rld) # Extract transformed OTU table
otu_table(fieldEco.vst) <- otu_table(fieldEco.vstMat, taxa_are_rows = TRUE) # Put transformed OTU table back into Phyloseq object
rm(fieldEco.dds,fieldEco.rld,fieldEco.vstMat) # Clean up
save(fieldEco.vst,file="intermediate_data/Phyloseq_fieldEco_vst.RData")
####### Register parallel backend #######
registerDoParallel(cores=4)

####### Figure S1: Ordination of ALL FIELD SAMPLES (leaves, roots, soil) #######
## Calculate weighted UniFrac distances
# ( replace negative values with 0 for distance calculations )
wUF.fieldEco.vst<-UniFrac(
  transform_sample_counts(fieldEco.vst,function(x) x<-ifelse(x<0,0,x)),
  weighted=TRUE,parallel=TRUE)
## Ordination: weighted UniFrac
PCoA.wUF.fieldEco.vst<-capscale(wUF.fieldEco.vst~1)
## Get inertia for top 3 PCoA axes: full dataset ##
PCoA.wUF.fieldEco.vst$CA$eig[1:3]/sum(PCoA.wUF.fieldEco.vst$CA$eig) 
## 50.0%, 12.8%, 5.0%
pdf(file="plots/Fig_S1_ordination_wUF_PCoA_vst_all_field_samples.pdf",width=9,height=9)
plot_ordination(fieldEco.vst,PCoA.wUF.fieldEco.vst,type='samples',axes=1:2,color='Site',shape='Type')+
  geom_point(size=3)+
  xlab("PCo1 [50.0%]")+ylab("PCo2 [12.8%]")+
  scale_color_manual(values=c(sitePalette[1],sitePalette[2],"mediumpurple","forest green",sitePalette[3]))+
  theme_classic()+ggtitle("Weighted UniFrac: all field samples")+
  theme(plot.title=element_text(size=28,face="bold"))+
  theme(axis.title.x=element_text(size=32,face="bold"),axis.text.x=element_text(size=26,face="bold"))+
  theme(axis.title.y=element_text(size=32,face="bold"),axis.text.y=element_text(size=26,face="bold"))+
  theme(legend.title= element_text(size=24),legend.text=element_text(size=22,face="bold"))+
  theme(legend.key.width=grid::unit(2,'lines'),legend.key.height=grid::unit(2,'lines'))
dev.off()
####### Table S1: Poor survival at Mil and Par #######
allsmd.leaf<-subset(as(sample_data(fieldEco.nobadOTUs.highcoverage.thresholded.CNC),"data.frame"),Genotype!='soil' & Type=='leaf'&Age!='Endogenous')
allsmd.root<-subset(as(sample_data(fieldEco.nobadOTUs.highcoverage.thresholded.CNC),"data.frame"),Genotype!='soil' & Type=='root'&Age!='Endogenous')
##Table S1: Poor survival at Mil and Par ##
sink("tables/Table_S1_poor_survival.txt")
print("Leaf:")
table(allsmd.leaf$Site,allsmd.leaf$Age)
print("Root:")
table(allsmd.root$Site,allsmd.root$Age)
sink()
rm(allsmd.leaf,allsmd.root) # clean up

####### Exclude samples from poor-survival sites; subset into leaves and roots for analysis #######
## Due to poor survival at gardens Par and Mil (resulting in unbalanced design), we exclude samples from those sites from further analysis: 
## Subset by organ type; keep endogenous plants for now
root3.withEndog.vst<-subset_samples(fieldEco.vst,Type=='root' & Site%in%c('Jam','Mah','Sil'))
leaf3.withEndog.vst<-subset_samples(fieldEco.vst,Type=='leaf' & Site%in%c('Jam','Mah','Sil'))
root3.withEndog<-subset_samples(fieldEco.nobadOTUs.highcoverage.thresholded.CNC,Type=='root' & Site%in%c('Jam','Mah','Sil')) %>%
  prune_taxa(taxa_sums(.)>0,.) # remove OTUs with 0 observations
leaf3.withEndog<-subset_samples(fieldEco.nobadOTUs.highcoverage.thresholded.CNC,Type=='leaf' & Site%in%c('Jam','Mah','Sil')) %>%
  prune_taxa(taxa_sums(.)>0,.) # remove OTUs with 0 observations
## These are the Phyloseq objects to be used in analyses going forward:
save(root3.withEndog.vst,file="intermediate_data/phylo_root3_withEndog_vst.RData")
save(root3.withEndog,file="intermediate_data/phylo_root3_withEndog.RData")
save(leaf3.withEndog.vst,file="intermediate_data/phylo_leaf3_withEndog_vst.RData")
save(leaf3.withEndog,file="intermediate_data/phylo_leaf3_withEndog.RData")
####### weighted UniFrac and PCoA: separately for leaf and root datasets #######
# replace negative values with 0s just for distance calculations
wUF.root3.withEndog.vst<-UniFrac(
  transform_sample_counts(root3.withEndog.vst,function(x) x<-ifelse(x<0,0,x)),
  weighted=TRUE,parallel=TRUE)
wUF.leaf3.withEndog.vst<-UniFrac(
  transform_sample_counts(leaf3.withEndog.vst,function(x) x<-ifelse(x<0,0,x)),
  weighted=TRUE,parallel=TRUE)

save(wUF.leaf3.withEndog.vst,file="intermediate_data/wUF_leaf3_wEndog_vst.RData")
save(wUF.root3.withEndog.vst,file="intermediate_data/wUF_root3_wEndog_vst.RData")

cap.wUF.root3.withEndog.vst<-capscale(wUF.root3.withEndog.vst~1,data=as(sample_data(root3.withEndog.vst),'data.frame'))
cap.wUF.leaf3.withEndog.vst<-capscale(wUF.leaf3.withEndog.vst~1,data=as(sample_data(leaf3.withEndog.vst),'data.frame'))

## Get inertia for top 3 PCoA axes: leaf ##
cap.wUF.leaf3.withEndog.vst$CA$eig[1:3]/sum(cap.wUF.leaf3.withEndog.vst$CA$eig) 
## 45.2%, 11.2%, 6.2%
cap.wUF.root3.withEndog.vst$CA$eig[1:3]/sum(cap.wUF.root3.withEndog.vst$CA$eig) 
## 29.6%, 16.1%, 10.6%

####### Fig. S5: Scree plots: weighted UniFrac #######
pdf(file="plots/Fig_S5a_screeplot_wUF_leaf.pdf")
barplot(cap.wUF.leaf3.withEndog.vst$CA$eig[1:20]/sum(cap.wUF.leaf3.withEndog.vst$CA$eig),
        main="Leaves: weighted UniFrac",xlab="Principal coordinate",ylab="Proportion Variance",
        cex.lab=2,cex.main=2.5,cex.axis=1.5,font.axis=2,font.main=2,font.lab=2,xaxt='n')
dev.off()
pdf(file="plots/Fig_S5b_screeplot_wUF_root.pdf")
barplot(cap.wUF.root3.withEndog.vst$CA$eig[1:20]/sum(cap.wUF.root3.withEndog.vst$CA$eig),
        main="Roots: weighted UniFrac",xlab="Principal coordinate",ylab="Proportion Variance",
        cex.lab=2,cex.main=2.5,cex.axis=1.5,font.axis=2,font.main=2,font.lab=2,xaxt='n')
dev.off()

####### How much variation is explained by the top 3 PCo axes? #######
sink("ordination_top3_cumulative_PVEs.txt")
print("cumulative percent variance explained by top 3 PCo:")
print("weighted UniFrac, roots:")
sum(cap.wUF.root3.withEndog.vst$CA$eig[1:3])/sum(cap.wUF.root3.withEndog.vst$CA$eig)  # 0.5635
print("weighted UniFrac, leaves:")
sum(cap.wUF.leaf3.withEndog.vst$CA$eig[1:3])/sum(cap.wUF.leaf3.withEndog.vst$CA$eig)  # 0.6266
print("Individual percent variance explained by top 3 PCo:")
print("weighted UniFrac, roots:")
(cap.wUF.root3.withEndog.vst$CA$eig[1:3])/sum(cap.wUF.root3.withEndog.vst$CA$eig) 
print("weighted UniFrac, leaves:")
(cap.wUF.leaf3.withEndog.vst$CA$eig[1:3])/sum(cap.wUF.leaf3.withEndog.vst$CA$eig) 
sink()

####### Fig. 4b: weighted UniFrac Ordination ~ Age #######
pdf(file="plots/Fig_4b_Ordination_Age_wUF1_2_leaf.pdf",width=9,height=9)
plot_ordination(leaf3.withEndog.vst,cap.wUF.leaf3.withEndog.vst,type="samples",axes=1:2,color="Age")+
  scale_colour_manual(values=agePalette)+
  geom_point(size=4,alpha=1)+
  xlab("PCo1 [45.2%]")+ylab("PCo2 [11.2%]")+
  ggtitle("Leaves")+theme_classic()+
  theme(plot.title = element_text(size=24, face="bold"))+
  theme(axis.title.x=element_text(size=24,face="bold"),axis.text.x=element_text(size=22,face="bold"))+
  theme(axis.title.y=element_text(size=24,face="bold"),axis.text.y=element_text(size=22,face="bold"))+
  theme(legend.title= element_text(size=24),legend.text=element_text(size=22))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

pdf(file="plots/Fig_4b_Ordination_Age_wUF2_3_root.pdf",width=9,height=9)
plot_ordination(root3.withEndog.vst,cap.wUF.root3.withEndog.vst,type="samples",axes=2:3,color="Age")+
  scale_colour_manual(values=agePalette,guide=FALSE)+
  geom_point(size=4,alpha=1)+
  xlab("PCo2 [16.1%]")+ylab("PCo3 [10.6%]")+
  ggtitle("Roots")+theme_classic()+
  theme(plot.title = element_text(size=24, face="bold"))+
  theme(axis.title.x=element_text(size=24,face="bold"),axis.text.x=element_text(size=22,face="bold"))+
  theme(axis.title.y=element_text(size=24,face="bold"),axis.text.y=element_text(size=22,face="bold"))+
  theme(legend.title= element_text(size=24),legend.text=element_text(size=22))+
  theme(legend.background = element_rect(fill="gray90", size=.5))+
  theme(axis.ticks=element_blank(),panel.grid.minor=element_blank())
dev.off()

####### Fig. 2b: weighted UniFrac Ordination ~ Site #######
pdf("plots/Fig_2b_Ordination_wUF1_2_Site_leaf.pdf",width=9,height=9)
plot_ordination(leaf3.withEndog.vst,cap.wUF.leaf3.withEndog.vst,type="samples",axes=1:2,color="Site")+
  scale_colour_manual(values=sitePalette)+
  geom_point(size=4,alpha=1)+
  xlab("PCo1 [45.3%]")+ylab("PCo2 [11.2%]")+
  ggtitle("Leaves")+theme_classic()+
  theme(plot.title = element_text(size=44, face="bold"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=40),legend.text=element_text(size=36,face="bold"))+
  theme(legend.key.height=unit(2.5,"lines"),legend.key.width=unit(2,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

pdf("plots/Fig_2b_Ordination_wUF1_2_Site_root.pdf",width=9,height=9)
plot_ordination(root3.withEndog.vst,cap.wUF.root3.withEndog.vst,type="samples",axes=1:2,color="Site")+
  scale_colour_manual(values=sitePalette,guide=FALSE)+
  geom_point(size=4,alpha=1)+
  xlab("PCo1 [29.6%]")+ylab("PCo2 [16.1%]")+
  ggtitle("Roots")+theme_classic()+
  theme(plot.title = element_text(size=44, face="bold"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))
dev.off()

####### Save major PCoA axes and Alpha diversity metrics for use in LMMs: #######
leaf3.smd.withEndog<-as(sample_data(leaf3.withEndog),'data.frame') %>%
  mutate(SampleID=row.names(.)) %>%
  merge(.,as.data.frame(cap.wUF.leaf3.withEndog.vst$CA$u[,1:3]),by.x='SampleID',by.y="row.names") %>%
  plyr::rename(replace=c('MDS1'='wUF1','MDS2'='wUF2','MDS3'='wUF3')) %>%
  merge(.,as.data.frame(estimate_richness(leaf3.withEndog,measures=c("Observed","Shannon","Chao1"))),by.x='SampleID',by.y='row.names') %>%
  mutate(expShannon=exp(Shannon)) # Shannon effective diversity

root3.smd.withEndog<-as(sample_data(root3.withEndog),'data.frame') %>%
  mutate(SampleID=row.names(.)) %>%
  merge(.,as.data.frame(cap.wUF.root3.withEndog.vst$CA$u[,1:3]),by.x='SampleID',by.y="row.names") %>%
  plyr::rename(replace=c('MDS1'='wUF1','MDS2'='wUF2','MDS3'='wUF3')) %>%
  merge(.,as.data.frame(estimate_richness(root3.withEndog,measures=c("Observed","Shannon","Chao1"))),by.x='SampleID',by.y='row.names') %>%
  mutate(expShannon=exp(Shannon)) # Shannon effective diversity

save(leaf3.smd.withEndog,file="intermediate_data/smd_leaf3_withEndog.RData")
save(root3.smd.withEndog,file="intermediate_data/smd_root3_withEndog.RData")

####### Filter out endogenous plants for analysis of experimental plants #######
leaf3<-subset_samples(leaf3.withEndog,Age!='endog') %>% 
  prune_taxa(taxa_sums(.)>0,.) # remove OTUs with 0 observations
root3<-subset_samples(root3.withEndog,Age!='endog') %>% 
  prune_taxa(taxa_sums(.)>0,.) # remove OTUs with 0 observations

ntaxa(leaf3) # 3611 of 3689 OTUs remain
ntaxa(root3) # 3677 of 3689 remain

save(leaf3,file="intermediate_data/phylo_leaf3.RData")
save(root3,file="intermediate_data/phylo_root3.RData")

####### How many OTUs unique to just leaves or roots? #######
leafonlyOTUs<-setdiff(taxa_names(leaf3),taxa_names(root3)) # 7 
rootonlyOTUs<-setdiff(taxa_names(root3),taxa_names(leaf3)) # 73

load("rel_abund/RelAbund_withEndog_otu.RData")
# how common are these?
filter(RelAbund.withEndog.otu, Taxon%in%leafonlyOTUs) %>% select(RA_leaf3) %>% sum # 6.36 e-6
filter(RelAbund.withEndog.otu, Taxon%in%rootonlyOTUs) %>% select(RA_root3) %>% sum # 0.00258

pdf(file="plots/Fig_S14a_leafonlyOTU_Abundance.pdf")
ggplot(RelAbund.withEndog.otu,aes(x=log(RA_leaf3)))+
  geom_density(aes(fill="all leaf OTUs\n",colour="all leaf OTUs\n"),alpha=0.3,position="identity")+
  geom_density(data=subset(RelAbund.withEndog.otu,Taxon%in%leafonlyOTUs),aes(x=log(RA_leaf3),colour="OTUs exclusively\nfound in leaves",fill="OTUs exclusively\nfound in leaves"),alpha=0.5,position="identity")+
  scale_x_continuous(labels=function(x) format(exp(x),digits=4),breaks=c(log(0.0000001),log(0.00001),log(0.001),log(0.1)))+
  scale_fill_manual(values=c("dark green","green"),name="")+
  scale_colour_manual(values=c("dark green","green"),name="")+
  theme_classic()+xlab("Relative abundance")+theme(legend.position=c(0.8,0.9))+
  theme(legend.title = element_text(size=20, face="bold"), legend.text = element_text(size=28,face="bold"))+
  theme(axis.title.x = element_text(size=36,face="bold"), axis.text.x = element_text(size=30,face="bold"))+
  theme(axis.title.y = element_text(size=36,face="bold"), axis.text.y = element_text(size=30,face="bold"))
dev.off()

pdf(file="plots/Fig_S14b_rootonlyOTU_Abundance.pdf")
ggplot(RelAbund.withEndog.otu,aes(x=log(RA_root3)))+
  geom_density(aes(fill="all root OTUs\n",colour="all root OTUs\n"),alpha=0.3,position="identity")+
  geom_density(data=subset(RelAbund.withEndog.otu,Taxon%in%rootonlyOTUs),aes(x=log(RA_root3),colour="OTUs exclusively\nfound in leaves",fill="OTUs exclusively\nfound in leaves"),alpha=0.5,position="identity")+
  scale_x_continuous(labels=function(x) format(exp(x),digits=4),breaks=c(log(0.0000001),log(0.00001),log(0.001),log(0.1)))+
  scale_fill_manual(values=c("dark grey","black"),name="")+
  scale_colour_manual(values=c("dark grey","black"),name="")+
  theme_classic()+xlab("Relative abundance")+theme(legend.position=c(0.8,0.9))+
  theme(legend.title = element_text(size=20, face="bold"), legend.text = element_text(size=28,face="bold"))+
  theme(axis.title.x = element_text(size=36,face="bold"), axis.text.x = element_text(size=30,face="bold"))+
  theme(axis.title.y = element_text(size=36,face="bold"), axis.text.y = element_text(size=30,face="bold"))
dev.off()

####### Conduct linear mixed models of alpha and beta diversity #######
# see LMMs.R

####### Fig. S1: OTU count distributions in experimental plants #######

# histograms of OTU counts
pdf(file="plots/Fig_S2_distributions_OTU_counts.pdf",width=9,height=6)
ggplot(rbind(data.frame("Observations"=taxa_sums(leaf3),"Organ"="leaf"),
             data.frame("Observations"=taxa_sums(root3),"Organ"="root")),
       aes(x=log(Observations),colour=Organ,fill=Organ))+
  scale_fill_manual(values=c("forest green","dark grey"))+
  scale_color_manual(values=c("forest green","dark grey"))+
  scale_x_continuous(labels=function(x) format(exp(x),digits=4),breaks=c(log(1),log(10),log(100),log(1000),log(10000),log(100000),log(1000000)))+
  ylab("Density")+xlab("Reads")+
  geom_density(alpha=0.3)+
  theme_classic()+
  theme(legend.title=element_blank(),legend.text=element_text(size=28,face="bold"),legend.position=c(0.9,0.8))+
  theme(axis.title.y=element_text(size=29,face="bold"),axis.text.y=element_text(size=26,face="bold"))+
  theme(axis.title.x=element_text(size=29,face="bold"),axis.text.x=element_text(size=26,face="bold"))
dev.off()

####### Fig. S6: Core microbiome analysis: Venn diagrams #######

## Core across sites:
pdf(file="plots/Fig_S6a_Venn_leaf3_Site.pdf")
grid.newpage()
venn.diagram(x=list(Jam=taxa_names(subset_samples(leaf3.withEndog,Site=='Jam')%>%prune_taxa(taxa_sums(.)>0,.)),
                    Mah=taxa_names(subset_samples(leaf3.withEndog,Site=='Mah')%>%prune_taxa(taxa_sums(.)>0,.)),
                    Sil=taxa_names(subset_samples(leaf3.withEndog,Site=='Sil')%>%prune_taxa(taxa_sums(.)>0,.))),
             cat.col=sitePalette[1:3],fill=sitePalette[1:3],alpha=c(0.6),cat.cex=2,cat.fontface="bold",cat.fontfamily="Helvetica",cat.dist=.08,
             label.col="white",cex=1.7,fontface="bold",fontfamily="Helvetica",
             main="leaf-associated OTUs",main.col="forest green",main.cex="2",main.fontface="bold",main.fontfamily="Helvetica",
             filename=NULL) %>% grid.draw()
dev.off()
pdf(file="plots/Fig_S6b_Venn_root3_Site.pdf")
grid.newpage()
venn.diagram(x=list(Jam=taxa_names(subset_samples(root3.withEndog,Site=='Jam')%>%prune_taxa(taxa_sums(.)>0,.)),
                    Mah=taxa_names(subset_samples(root3.withEndog,Site=='Mah')%>%prune_taxa(taxa_sums(.)>0,.)),
                    Sil=taxa_names(subset_samples(root3.withEndog,Site=='Sil')%>%prune_taxa(taxa_sums(.)>0,.))),
             cat.col=sitePalette[1:3],fill=sitePalette[1:3],alpha=c(0.6),cat.cex=2,cat.fontface="bold",cat.fontfamily="Helvetica",cat.dist=.08,
             label.col="white",cex=1.7,fontface="bold",fontfamily="Helvetica",
             main="root-associated OTUs",main.col="grey",main.cex="2",main.fontface="bold",main.fontfamily="Helvetica",
             filename=NULL) %>% grid.draw()
dev.off()

## Core across Genotypes:
pdf(file="plots/Fig_S6b_Venn_leaf3_Geno.pdf")
grid.newpage()
venn.diagram(x=list(JAM=taxa_names(subset_samples(leaf3.withEndog,Genotype=='JAM')%>%prune_taxa(taxa_sums(.)>0,.)),
                    MAH=taxa_names(subset_samples(leaf3.withEndog,Genotype=='MAH')%>%prune_taxa(taxa_sums(.)>0,.)),
                    MIL=taxa_names(subset_samples(leaf3.withEndog,Genotype=='MIL')%>%prune_taxa(taxa_sums(.)>0,.)),
                    PAR=taxa_names(subset_samples(leaf3.withEndog,Genotype=='PAR')%>%prune_taxa(taxa_sums(.)>0,.)),
                    SIL=taxa_names(subset_samples(leaf3.withEndog,Genotype=='SIL')%>%prune_taxa(taxa_sums(.)>0,.))),
             cat.col=popPalette[1:5],fill=popPalette[1:5],alpha=c(0.6),cat.cex=2,cat.fontface="bold",cat.fontfamily="Helvetica",cat.dist=.22,
             label.col="white",cex=1.7,fontface="bold",fontfamily="Helvetica",main.just=c(2,5),
             main="leaf-associated OTUs \nby host genotype\n",main.col="forest green",main.cex="2",main.fontface="bold",main.fontfamily="Helvetica",
             filename=NULL) %>% grid.draw()
dev.off()
pdf(file="plots/Fig_S6b_Venn_root3_Geno.pdf")
grid.newpage()
venn.diagram(x=list(JAM=taxa_names(subset_samples(root3.withEndog,Genotype=='JAM')%>%prune_taxa(taxa_sums(.)>0,.)),
                    MAH=taxa_names(subset_samples(root3.withEndog,Genotype=='MAH')%>%prune_taxa(taxa_sums(.)>0,.)),
                    MIL=taxa_names(subset_samples(root3.withEndog,Genotype=='MIL')%>%prune_taxa(taxa_sums(.)>0,.)),
                    PAR=taxa_names(subset_samples(root3.withEndog,Genotype=='PAR')%>%prune_taxa(taxa_sums(.)>0,.)),
                    SIL=taxa_names(subset_samples(root3.withEndog,Genotype=='SIL')%>%prune_taxa(taxa_sums(.)>0,.))),
             cat.col=popPalette[1:5],fill=popPalette[1:5],alpha=c(0.6),cat.cex=2,cat.fontface="bold",cat.fontfamily="Helvetica",cat.dist=.22,
             label.col="white",cex=1.7,fontface="bold",fontfamily="Helvetica",main.just=c(2,5),
             main="root-associated OTUs \nby host genotype\n",main.col="grey",main.cex="2",main.fontface="bold",main.fontfamily="Helvetica",
             filename=NULL) %>% grid.draw()
dev.off()

## Core across ages:
pdf(file="plots/Fig_S6b_Venn_leaf3_Age.pdf")
grid.newpage()
venn.diagram(x=list(age2=taxa_names(subset_samples(leaf3.withEndog,Age=='2')%>%prune_taxa(taxa_sums(.)>0,.)),
                    age3=taxa_names(subset_samples(leaf3.withEndog,Age=='3')%>%prune_taxa(taxa_sums(.)>0,.)),
                    age4=taxa_names(subset_samples(leaf3.withEndog,Age=='4')%>%prune_taxa(taxa_sums(.)>0,.))),
             cat.col=agePalette[1:3],fill=agePalette[1:3],alpha=c(0.5),cat.cex=2,cat.fontface="bold",cat.fontfamily="Helvetica",cat.dist=.08,
             label.col="white",cex=1.7,fontface="bold",fontfamily="Helvetica",
             main="leaf-associated OTUs by host Age",main.col="forest green",main.cex="2",main.fontface="bold",main.fontfamily="Helvetica",
             filename=NULL) %>% grid.draw()
dev.off()
pdf(file="plots/Fig_S6b_Venn_root3_Age.pdf")
grid.newpage()
venn.diagram(x=list(age2=taxa_names(subset_samples(root3.withEndog,Age=='2')%>%prune_taxa(taxa_sums(.)>0,.)),
                    age3=taxa_names(subset_samples(root3.withEndog,Age=='3')%>%prune_taxa(taxa_sums(.)>0,.)),
                    age4=taxa_names(subset_samples(root3.withEndog,Age=='4')%>%prune_taxa(taxa_sums(.)>0,.))),
             cat.col=agePalette[1:3],fill=agePalette[1:3],alpha=c(0.5),cat.cex=2,cat.fontface="bold",cat.fontfamily="Helvetica",cat.dist=.08,
             label.col="white",cex=1.7,fontface="bold",fontfamily="Helvetica",
             main="root-associated OTUs by host Age",main.col="grey",main.cex="2",main.fontface="bold",main.fontfamily="Helvetica",
             filename=NULL) %>% grid.draw()
dev.off()

####### Fig. S2: Soil nutrition data #######
soil<-as.data.frame(read.delim("raw_data//soildata.txt",sep='\t',header=TRUE))
soil<-filter(soil,!(Site%in%c('Ald','Par'))) # focus on the main 3 sites: Jam, Mah, Sil
soil$Site<-factor(soil$Site)

soilPCA<-capscale(select(soil, -Site)~1,distance="euclidean")
summary(soilPCA) # PC1 = 69.97 % variation; # PC2 = 27.52 % variation
soil<-merge(soil,scores(soilPCA,display="sites",choices=c(1:3)),by='row.names') # save first 3 PCo
pdf(file="plots/Fig_S3b_soilPCA.pdf",width=9,height=9)
  ggplot(soil,aes(x=MDS1,y=MDS2,color=Site))+
    geom_point(size=5)+
    scale_color_manual(values=c(sitePalette[1:3]))+
    xlab("PC1 [70.0%]")+ylab("PC2 [27.5%]")+
    theme_classic()+theme(plot.title = element_text(size=24, face="bold"))+
    theme(axis.title.x=element_text(size=44,face="bold"),axis.text.x=element_text(size=36,face="bold"))+
    theme(axis.title.y=element_text(size=44,face="bold"),axis.text.y=element_text(size=36,face="bold"))+
    theme(legend.title= element_text(size=44),legend.text=element_text(size=36))+
    theme(legend.key.height=unit(2.5,'lines'),legend.key.width=unit(2,'lines'))+
    theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

soil.melted<-melt(soil,id.vars=c("Row.names","Site"))
pdf(file="plots/Fig_S3a_Soil_boxplots.pdf",width=9,height=9)
update_geom_defaults("point", list(colour = NULL))
filter(soil.melted,Site!='Par',!(variable%in%c('MDS1','MDS2','MDS3')))%>%
  ggplot(.,aes(x=Site,y=value,color=Site))+
    facet_wrap(~variable,scales="free")+
    scale_color_manual(values=sitePalette,name="Site:   ")+
    geom_boxplot(size=1.5)+
    theme_classic()+theme(legend.direction="horizontal")+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_text(size=22,face="bold"))+
    theme(legend.title= element_text(size=24),legend.text=element_text(size=22,face="bold"))+
    theme(legend.background = element_rect(fill="gray90", size=.5))+
    theme(legend.position='top',legend.key.height=unit(2,'lines'),legend.key.width=unit(2,'lines'))+
    theme(strip.background=element_rect(fill="gray90"),strip.text=element_text(size=20,face="bold"))
update_geom_defaults("point", list(colour = "black"))
dev.off()

####### Fig. 1a: map #######
# make inset
pdf(file="plots/Fig_1a_inset.pdf",width=7,height=7)
map(database="state",regions=c("Montana","Idaho", "Wyoming","Oregon","Washington"),
    xlim=c(-117.25,-111.05),ylim=c(42,49),fill=TRUE,col="white")
map(database="state",regions=c("Idaho"),add=TRUE,fill=TRUE,col="gray90")
polygon(x=c(-113,-113,-115,-115),y=c(43.7,46,46,43.7),border="black",lwd=2)
dev.off()

# load site coordinates
coords<-as.data.frame(read.table("raw_data/site_coords_Ecotypes.txt",sep="\t",header=TRUE))

# make regional map
pdf(file="plots/Fig_1a.pdf",width=5,height=7)
map(database="county",regions=c("Idaho"),xlim=c(-115,-113),ylim=c(43.7,46))
polygon(x=c(-113,-113,-115,-115),y=c(43.7,46,46,43.7),border="black",lwd=8)
#points(x=subset(coords,Site%in%c("Jam","Mah","Sil"))$Lon, y=subset(coords,Site%in%c("Jam","Mah","Sil"))$Lat, pch=17,cex=3)
points(x=subset(coords,Site=="Jam")$Lon, y=subset(coords,Site=="Jam")$Lat, pch=17,cex=3,col=sitePalette[1])
points(x=subset(coords,Site=="Mah")$Lon, y=subset(coords,Site=="Mah")$Lat, pch=17,cex=3,col=sitePalette[2])
points(x=subset(coords,Site=="Sil")$Lon, y=subset(coords,Site=="Sil")$Lat, pch=17,cex=3,col=sitePalette[3])
points(x=subset(coords,Site%in%c("Mil","Par"))$Lon, y=subset(coords,Site%in%c("Mil","Par"))$Lat, pch=19,cex=3,col="dark grey")
text(x=c(-113.8,-113.8,-114.6,-113.3,-114.5),
     y=c(45.1,44.05,45.02,44.5,44.5),
     labels=c("Jam","Mah","Sil","Mil","Par"),
    col=c(sitePalette,"dark grey","dark grey"),cex=2,font=2,family="Helvetica")
# add scale bar (50 km calculated using http://www.nhc.noaa.gov/gccalc.shtml)
polygon(x=c(-113.84,-113.2),y=c(45.5,45.5),lwd=8,border="black")
# text(x=-113.5,y=45.6,labels="50 km",cex=1.3,font=2)
dev.off()

####### Estimate fold changes of individual taxa #######
# see ./foldchange/foldchange.R

####### Root-leaf overlap: 99%-binned OTU analyses  #######
# see otu99pct/otu99pct.R

####### Parallel analysis using UNweighted UniFrac metrics #######
# see uUF/uUF.R

####### Analysis of glucosinolate profiles #######
# see GLS/GLS.R

####### Figure 1e: phylum barplots #######
load('higher_tax_levels/phyloseq_fullEco_CNC_unfiltered_phy.RData')

arrange(RelAbund.withEndog.phy,-RA_root5) %>% 
  mutate(Taxon=ordered(Taxon,levels=Taxon)) %>%
  select(Taxon) %>%
  slice(1:12) %>% 
  print
  
fieldEco.phy<-subset_samples(fullEco.phy,Site%in%c('Jam','Mah','Sil') & Type!='soil') %>%
  prune_taxa(taxa_sums(.)>0,.) # remove taxa with 0 observations

# Find top 12 phyla by abundance
sort(taxa_sums(fieldEco.phy),decreasing=TRUE)[1:12] %>% names() -> top12phyla
setdiff(taxa_names(fieldEco.phy),top12phyla) -> rarephyla

# Consolidate rare phyla into "Low abundance" category
merged.fieldEco.phy<-merge_taxa(fieldEco.phy,rarephyla,archetype="TM7") %>%
  prune_taxa(taxa_sums(.)>0,.)
ntaxa(merged.fieldEco.phy) # 13 taxa remain, as expected
taxa_names(merged.fieldEco.phy)

# Extract phylum counts & sample metadata into a data frame
phydata<-as(otu_table(merged.fieldEco.phy),'matrix') %>% t() %>% as.data.frame() %>% 
  mutate('SampleID'=factor(row.names(.))) %>% 
  merge(.,as(sample_data(merged.fieldEco.phy),'data.frame')%>%select(SampleID,Type,Site),by='SampleID') %>% # merge with sample metadata
  gather(key=Phylum,value=Abundance,2:14) %>% # melt data frame
  group_by(SampleID) %>% mutate(RelAbund=Abundance/sum(Abundance)) %>%   # transform into relative abundance in each sample:
  ungroup %>% as.data.frame %>%
  mutate(Phylum=ifelse(Phylum=='TM7','Low abundance',Phylum)) # rename TM7 as 'Low abundance' (TM7 is a placeholder for merged rare phyla)


phyPalette<-c("#000000",
             "#B15928",
             "#33A02C",
             "#B2DF8A",
             "#1F78B4",
             "#FB9A99",
             "#FDBF6F",
             "#E31A1C",
             "#6A3D9A",
             "#CAB2D6",
             "#FF7F00",
             "#FFFF99",
             "#A6CEE3")

# Plot for leaves:
phyleaf<-
  filter(phydata,Type=='leaf') %>%
  mutate(SampleID=factor(SampleID)) %>%
  mutate(Phylum=ordered(Phylum,levels=rev(c('Low abundance',top12phyla[12:1])))) %>%
  arrange(Phylum) %>%
ggplot(.,aes(x=SampleID,y=100*RelAbund,fill=Phylum,color=Phylum)) +
  geom_bar(stat='identity') +
  facet_wrap(~Site,scales='free_x',switch='x')+
  theme_classic()+ 
  ylab("Relative abundance (%)")+xlab("Site")+
  ggtitle("Leaves")+
  scale_fill_manual(values=phyPalette[13:1],guide=FALSE)+
  scale_color_manual(values=phyPalette[13:1],guide=FALSE)+
  scale_x_discrete(breaks=NULL)+
  theme(plot.title=element_text(size=32,face='bold'))+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=26,face='bold',vjust=1.5))+
  theme(axis.title.x=element_text(size=30,face='bold'))+
  theme(axis.title.y=element_text(size=30,face='bold'),axis.text.y=element_text(size=24,face='bold'))

# Plot for roots:
phyroot<-
  filter(phydata,Type=='root') %>%
  mutate(SampleID=factor(SampleID)) %>%
  mutate(Phylum=ordered(Phylum,levels=rev(c('Low abundance',top12phyla[12:1])))) %>%
  arrange(Phylum) %>%
  ggplot(.,aes(x=SampleID,y=100*RelAbund,fill=Phylum,color=Phylum)) +
  geom_bar(stat='identity') +
  facet_wrap(~Site,scales='free_x',switch='x')+
  theme_classic()+ 
  ylab("Relative abundance (%)")+xlab("Site")+
  ggtitle("Roots")+
  scale_fill_manual(values=phyPalette[13:1],guide=FALSE)+
  scale_color_manual(values=phyPalette[13:1],guide=FALSE)+
  scale_x_discrete(breaks=NULL)+
  scale_y_continuous(breaks=NULL)+
  theme(plot.title=element_text(size=32,face='bold'))+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=26,face='bold',vjust=1.5))+
  theme(axis.title.x=element_text(size=30,face='bold'))+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank())

# Plot for legend only (combine with data figures in Inkscape) 
pdf(file="plots/Fig_1e_legend.pdf",width=5,height=9)
mutate(phydata,Phylum=ordered(Phylum,levels=rev(c('Low abundance',top12phyla[12:1])))) %>%
  arrange(Phylum) %>%
ggplot(.,aes(x=Type,y=RelAbund,fill=Phylum))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=phyPalette[13:1])+
  guides(fill=guide_legend(reverse=TRUE))+
  theme(legend.title=element_text(size=20,face='bold'),legend.text=element_text(size=14,face='bold'))
dev.off()

pdf("plots/Fig_1e_leaf_phyla_barplot.pdf",width=9,height=9)
phyleaf
dev.off()

pdf("plots/Fig_1e_root_phyla_barplot.pdf",width=9,height=9)
phyroot
dev.off()

####### Save image #######
save.image(paste0("main_code_image_",date(),".RData"))
