### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
#####  Consolidate count data at higher taxonomic levels instead of OTU level

####### Clear workspace ########
rm(list=ls())

####### Load source file #######
source('ecotypes_source.R')

####### Define zero-tolerant function for calculating geometric means (credit: P.J. McMurdie) #######
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

####### Load data file (has NOT been thresholded or copy number corrected) ######
load("intermediate_data/fullEco_preThreshold_preCNC.RData")

####### Copy number correction at OTU level (pre-filtering) #######
## Load copy number estimates:
CNCest<-as.data.frame(read.table('raw_data/Eco_Field_copynumest_forR.txt',sep='\t',header=TRUE))

fullOTU.unfiltered.CNC<-as.data.frame(as(otu_table(fullEco.nobadOTUs.highcoverage),'matrix')) # copy OTU table from thresholded Phyloseq object
fullOTU.unfiltered.CNC<-merge(CNCest,fullOTU.unfiltered.CNC,by.y='row.names',by.x='OTU_ID') # merge copy number estimates with OTU table
rownames(fullOTU.unfiltered.CNC)<-fullOTU.unfiltered.CNC$OTU_ID # restore OTU_IDs as row names
fullOTU.unfiltered.CNC<-fullOTU.unfiltered.CNC[,3:dim(fullOTU.unfiltered.CNC)[2]]/fullOTU.unfiltered.CNC[,2] # divide OTU counts by corresponding copy number estimate (Col 1 = OTU_ID, Col 2 = copy number estimates)
fullOTU.unfiltered.CNC<-data.matrix(ceiling(fullOTU.unfiltered.CNC)) # round OTU counts up to nearest integer

## Replace OTU table of master Phyloseq objects with copy number corrected versions:
# do for fullEco dataset:
fullEco.nobadOTUs.highcoverage.CNC<-fullEco.nobadOTUs.highcoverage # copy Phyloseq object
otu_table(fullEco.nobadOTUs.highcoverage.CNC)<-otu_table(fullOTU.unfiltered.CNC,taxa_are_rows=TRUE) # replace OTU table with copy number-corrected version
save(fullEco.nobadOTUs.highcoverage.CNC,file="higher_tax_levels/phyloseq_fullEco_CNC_unfiltered.RData")

####### Tax_glom() to pool OTUs by higher taxonomic levels #######
fullEco.fam<-tax_glom(fullEco.nobadOTUs.highcoverage.CNC,taxrank="Family")
taxa_names(fullEco.fam)<-as.data.frame(tax_table(fullEco.fam))$Family # update taxa names
save(fullEco.fam,file="higher_tax_levels/phyloseq_fullEco_CNC_unfiltered_fam.RData")

fullEco.ord<-tax_glom(fullEco.nobadOTUs.highcoverage.CNC,taxrank="Order")
taxa_names(fullEco.ord)<-as.data.frame(tax_table(fullEco.ord))$Order # update taxa names
save(fullEco.ord,file="higher_tax_levels/phyloseq_fullEco_CNC_unfiltered_ord.RData")

fullEco.cla<-tax_glom(fullEco.nobadOTUs.highcoverage.CNC,taxrank="Class")
taxa_names(fullEco.cla)<-as.data.frame(tax_table(fullEco.cla))$Class # update taxa names
save(fullEco.cla,file="higher_tax_levels/phyloseq_fullEco_CNC_unfiltered_cla.RData")

fullEco.phy<-tax_glom(fullEco.nobadOTUs.highcoverage.CNC,taxrank="Phylum")
taxa_names(fullEco.phy)<-as.data.frame(tax_table(fullEco.phy))$Phylum # update taxa names
save(fullEco.phy,file="higher_tax_levels/phyloseq_fullEco_CNC_unfiltered_phy.RData")

# To load these files again after they were first created (if needed):
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_fam.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_ord.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_cla.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_phy.RData")

####### variance stabilizing transformation at each taxonomic level-- for heritability estimates using LMMs #######
## Apply variance stabilizing transformation to *all FIELD samples* (same as for OTUs in main_code)
fieldEco.phy.vst<-subset_samples(fullEco.phy, Site!='Duke') %>% prune_taxa(taxa_sums(.)>0,.) # copy Phyloseq object while removing taxa with no observations
fieldEco.phy.dds<-phyloseq_to_deseq2(fieldEco.phy.vst,~Type+Site+Type*Site) # Change to DESeq2 object
fieldEco.phy.dds = estimateSizeFactors(fieldEco.phy.dds, geoMeans = apply(counts(fieldEco.phy.dds), 1, gm_mean)) # Estimate size factors using geometric means 
fieldEco.phy.rld <- DESeq2::varianceStabilizingTransformation(fieldEco.phy.dds, blind = FALSE) # Apply variance stabilizing transformation
fieldEco.phy.vstMat <- GenomicRanges::assay(fieldEco.phy.rld) # Extract transformed OTU table
otu_table(fieldEco.phy.vst) <- otu_table(fieldEco.phy.vstMat, taxa_are_rows = TRUE) # Put transformed OTU table back into Phyloseq object
rm(fieldEco.phy.dds,fieldEco.phy.rld,fieldEco.phy.vstMat) # Clean up
print("Phylum vst transformation complete")
save(fieldEco.phy.vst,file="higher_tax_levels/vst_fieldEco_phy_CNC_unfiltered.RData")

## Apply variance stabilizing transformation to *all FIELD samples* (same as for OTUs in main_code)
fieldEco.cla.vst<-subset_samples(fullEco.cla, Site!='Duke') %>% prune_taxa(taxa_sums(.)>0,.) # copy Phyloseq object while removing taxa with no observations
fieldEco.cla.dds<-phyloseq_to_deseq2(fieldEco.cla.vst,~Type+Site+Type*Site) # Change to DESeq2 object
fieldEco.cla.dds = estimateSizeFactors(fieldEco.cla.dds, geoMeans = apply(counts(fieldEco.cla.dds), 1, gm_mean)) # Estimate size factors using geometric means 
fieldEco.cla.rld <- DESeq2::varianceStabilizingTransformation(fieldEco.cla.dds, blind = FALSE) # Apply variance stabilizing transformation
fieldEco.cla.vstMat <- GenomicRanges::assay(fieldEco.cla.rld) # Extract transformed OTU table
otu_table(fieldEco.cla.vst) <- otu_table(fieldEco.cla.vstMat, taxa_are_rows = TRUE) # Put transformed OTU table back into Phyloseq object
rm(fieldEco.cla.dds,fieldEco.cla.rld,fieldEco.cla.vstMat) # Clean up
print("Class vst transformation complete")
save(fieldEco.cla.vst,file="higher_tax_levels/vst_fieldEco_cla_CNC_unfiltered.RData")

## Apply variance stabilizing transformation to *all FIELD samples* (same as for OTUs in main_code)
fieldEco.ord.vst<-subset_samples(fullEco.ord, Site!='Duke') %>% prune_taxa(taxa_sums(.)>0,.) # copy Phyloseq object while removing taxa with no observations
fieldEco.ord.dds<-phyloseq_to_deseq2(fieldEco.ord.vst,~Type+Site+Type*Site) # Change to DESeq2 object
fieldEco.ord.dds = estimateSizeFactors(fieldEco.ord.dds, geoMeans = apply(counts(fieldEco.ord.dds), 1, gm_mean)) # Estimate size factors using geometric means 
fieldEco.ord.rld <- DESeq2::varianceStabilizingTransformation(fieldEco.ord.dds, blind = FALSE) # Apply variance stabilizing transformation
fieldEco.ord.vstMat <- GenomicRanges::assay(fieldEco.ord.rld) # Extract transformed OTU table
otu_table(fieldEco.ord.vst) <- otu_table(fieldEco.ord.vstMat, taxa_are_rows = TRUE) # Put transformed OTU table back into Phyloseq object
rm(fieldEco.ord.dds,fieldEco.ord.rld,fieldEco.ord.vstMat) # Clean up
print("Order vst transformation complete")
save(fieldEco.ord.vst,file="higher_tax_levels/vst_fieldEco_ord_CNC_unfiltered.RData")

## Apply variance stabilizing transformation to *all FIELD samples* (same as for OTUs in main_code)
fieldEco.fam.vst<-subset_samples(fullEco.fam, Site!='Duke') %>% prune_taxa(taxa_sums(.)>0,.) # copy Phyloseq object while removing taxa with no observations
fieldEco.fam.dds<-phyloseq_to_deseq2(fieldEco.fam.vst,~Type+Site+Type*Site) # Change to DESeq2 object
fieldEco.fam.dds = estimateSizeFactors(fieldEco.fam.dds, geoMeans = apply(counts(fieldEco.fam.dds), 1, gm_mean)) # Estimate size factors using geometric means 
fieldEco.fam.rld <- DESeq2::varianceStabilizingTransformation(fieldEco.fam.dds, blind = FALSE) # Apply variance stabilizing transformation
fieldEco.fam.vstMat <- GenomicRanges::assay(fieldEco.fam.rld) # Extract transformed OTU table
otu_table(fieldEco.fam.vst) <- otu_table(fieldEco.fam.vstMat, taxa_are_rows = TRUE) # Put transformed OTU table back into Phyloseq object
rm(fieldEco.fam.dds,fieldEco.fam.rld,fieldEco.fam.vstMat) # Clean up
print("Family vst transformation complete")
save(fieldEco.fam.vst,file="higher_tax_levels/vst_fieldEco_fam_CNC_unfiltered.RData")

####### Save image #######
save.image(paste0("higher_tax_levels/image_",date()))