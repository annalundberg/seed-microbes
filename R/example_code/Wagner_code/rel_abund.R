### Analyses by Maggie R. Wagner
### maggie.r.wagner@gmail.com
#####  Calculate relative abundances of OTUs & higher taxonomic levels in various subsets of samples 

####### Clear workspace ########
rm(list=ls())

####### Load source file #######
source('ecotypes_source.R')

####### Load data files: OTU level ######

load("intermediate_data/fieldEco_cleaned_CNC.RData")
fieldEco<-fieldEco.nobadOTUs.highcoverage.thresholded.CNC; rm(fieldEco.nobadOTUs.highcoverage.thresholded.CNC) # rename for convenience
load("intermediate_data/fullEco_cleaned_CNC.RData")
fullEco<-fullEco.nobadOTUs.highcoverage.thresholded.CNC; rm(fullEco.nobadOTUs.highcoverage.thresholded.CNC) # rename for convenience

####### Load data files: higher tax. levels #######
## Note: unlike the OTU-level data, these were NOT pre-filtered/thresholded for low abundance. But they WERE copy-number corrected
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_fam.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_ord.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_cla.RData")
load("higher_tax_levels//phyloseq_fullEco_CNC_unfiltered_phy.RData")

####### OTU level: Make dataframe of relative OTU abundances in various subsets of samples -- with Endogenous plants #######
RelAbund.withEndog.otu<-data.frame("Taxon"=taxa_names(fullEco),"Level"='otu') %>%
  mutate("RA_leaf5"=taxa_sums(subset_samples(fullEco, Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco, Type=='leaf')))) %>%
  mutate("RA_root5"=taxa_sums(subset_samples(fullEco, Type=='root'))/sum(taxa_sums(subset_samples(fullEco, Type=='root')))) %>%
  mutate("RA_leaf3"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &  Type=='leaf')))) %>%
  mutate("RA_root3"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &  Type=='root')))) %>%
  mutate("RA_leaf_Jam"=taxa_sums(subset_samples(fullEco,Site=='Jam' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site=='Jam' &  Type=='leaf')))) %>%
  mutate("RA_root_Jam"=taxa_sums(subset_samples(fullEco,Site=='Jam' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site=='Jam' &  Type=='root')))) %>%
  mutate("RA_leaf_Mah"=taxa_sums(subset_samples(fullEco,Site=='Mah' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site=='Mah' &  Type=='leaf')))) %>%
  mutate("RA_root_Mah"=taxa_sums(subset_samples(fullEco,Site=='Mah' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site=='Mah' &  Type=='root')))) %>%
  mutate("RA_leaf_Sil"=taxa_sums(subset_samples(fullEco,Site=='Sil' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site=='Sil' &  Type=='leaf')))) %>%
  mutate("RA_root_Sil"=taxa_sums(subset_samples(fullEco,Site=='Sil' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site=='Sil' &  Type=='root')))) %>%
  mutate("RA_root_Par"=taxa_sums(subset_samples(fullEco,Site=='Par' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site=='Par' & Type=='root')))) %>%
  mutate("RA_root_Mil"=taxa_sums(subset_samples(fullEco,Site=='Mil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site=='Mil' & Type=='root')))) %>%
  mutate("RA_leaf_Par"=taxa_sums(subset_samples(fullEco,Site=='Par' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site=='Par' & Type=='leaf')))) %>%
  mutate("RA_leaf_Mil"=taxa_sums(subset_samples(fullEco,Site=='Mil' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site=='Mil' & Type=='leaf')))) %>%
  mutate("RA_root_GRHsoil"=taxa_sums(subset_samples(fullEco,Treatment=='GRHsoil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Treatment=='GRHsoil' & Type=='root')))) %>%
  mutate("RA_soil_Jam"=taxa_sums(subset_samples(fullEco,Site=='Jam' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco,Site=='Jam' & Type=='soil')))) %>%
  mutate("RA_soil_Mah"=taxa_sums(subset_samples(fullEco,Site=='Mah' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco,Site=='Mah' & Type=='soil')))) %>%
  mutate("RA_soil_Sil"=taxa_sums(subset_samples(fullEco,Site=='Sil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco,Site=='Sil' & Type=='soil')))) %>%
  mutate("RA_soil_Par"=taxa_sums(subset_samples(fullEco,Site=='Par' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco,Site=='Par' & Type=='soil')))) %>%
  mutate("RA_soil_Mil"=taxa_sums(subset_samples(fullEco,Site=='Mil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco,Site=='Mil' & Type=='soil')))) %>%
  mutate("RA_soil_GRH"=taxa_sums(subset_samples(fullEco,Treatment=='GRHsoil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco,Treatment=='GRHsoil' & Type=='soil')))) %>%
  mutate("RA_leaf3_age2"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') & Age=='2' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age2"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') & Age=='2'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age3"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') & Age=='3' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age3"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') & Age=='3'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age4"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') & Age=='4' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age4"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') & Age=='4'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='root')))) %>%
  mutate("RA_leaf3_endog"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='leaf')))) %>%
  mutate("RA_root3_endog"=taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='root'))))
save(RelAbund.withEndog.otu,file="rel_abund/RelAbund_withEndog_otu.RData")
####### Family level: Make dataframe of relative Family abundances in various subsets of samples -- with Endogenous plants #######
RelAbund.withEndog.fam<-data.frame("Taxon"=taxa_names(fullEco.fam),"Level"='fam') %>%
  mutate("RA_leaf5"=taxa_sums(subset_samples(fullEco.fam, Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam, Type=='leaf')))) %>%
  mutate("RA_root5"=taxa_sums(subset_samples(fullEco.fam, Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam, Type=='root')))) %>%
  mutate("RA_leaf3"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &  Type=='leaf')))) %>%
  mutate("RA_root3"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &  Type=='root')))) %>%
  mutate("RA_leaf_Jam"=taxa_sums(subset_samples(fullEco.fam,Site=='Jam' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Jam' &  Type=='leaf')))) %>%
  mutate("RA_root_Jam"=taxa_sums(subset_samples(fullEco.fam,Site=='Jam' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Jam' &  Type=='root')))) %>%
  mutate("RA_leaf_Mah"=taxa_sums(subset_samples(fullEco.fam,Site=='Mah' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Mah' &  Type=='leaf')))) %>%
  mutate("RA_root_Mah"=taxa_sums(subset_samples(fullEco.fam,Site=='Mah' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Mah' &  Type=='root')))) %>%
  mutate("RA_leaf_Sil"=taxa_sums(subset_samples(fullEco.fam,Site=='Sil' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Sil' &  Type=='leaf')))) %>%
  mutate("RA_root_Sil"=taxa_sums(subset_samples(fullEco.fam,Site=='Sil' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Sil' &  Type=='root')))) %>%
  mutate("RA_root_Par"=taxa_sums(subset_samples(fullEco.fam,Site=='Par' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Par' & Type=='root')))) %>%
  mutate("RA_root_Mil"=taxa_sums(subset_samples(fullEco.fam,Site=='Mil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Mil' & Type=='root')))) %>%
  mutate("RA_leaf_Par"=taxa_sums(subset_samples(fullEco.fam,Site=='Par' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Par' & Type=='leaf')))) %>%
  mutate("RA_leaf_Mil"=taxa_sums(subset_samples(fullEco.fam,Site=='Mil' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Mil' & Type=='leaf')))) %>%
  mutate("RA_root_GRHsoil"=taxa_sums(subset_samples(fullEco.fam,Treatment=='GRHsoil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Treatment=='GRHsoil' & Type=='root')))) %>%
  mutate("RA_soil_Jam"=taxa_sums(subset_samples(fullEco.fam,Site=='Jam' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Jam' & Type=='soil')))) %>%
  mutate("RA_soil_Mah"=taxa_sums(subset_samples(fullEco.fam,Site=='Mah' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Mah' & Type=='soil')))) %>%
  mutate("RA_soil_Sil"=taxa_sums(subset_samples(fullEco.fam,Site=='Sil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Sil' & Type=='soil')))) %>%
  mutate("RA_soil_Par"=taxa_sums(subset_samples(fullEco.fam,Site=='Par' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Par' & Type=='soil')))) %>%
  mutate("RA_soil_Mil"=taxa_sums(subset_samples(fullEco.fam,Site=='Mil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.fam,Site=='Mil' & Type=='soil')))) %>%
  mutate("RA_soil_GRH"=taxa_sums(subset_samples(fullEco.fam,Treatment=='GRHsoil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.fam,Treatment=='GRHsoil' & Type=='soil')))) %>%
  mutate("RA_leaf3_age2"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age=='2' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age2"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age=='2'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age3"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age=='3' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age3"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age=='3'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age4"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age=='4' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age4"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age=='4'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='root')))) %>%
  mutate("RA_leaf3_endog"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='leaf')))) %>%
  mutate("RA_root3_endog"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='root')))) %>%
  mutate("RA_root3_exp"=taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') & Age!='endog' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site %in% c('Jam','Mah','Sil') &Age!='endog'  & Type=='root')))) %>%
  mutate("RA_root5_endog"=taxa_sums(subset_samples(fullEco.fam,Site!='Duke' & Age=='endog' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site!='Duke' & Age=='endog'  & Type=='root')))) %>%
  mutate("RA_root5_exp"=taxa_sums(subset_samples(fullEco.fam,Site!='Duke' & Age!='endog' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.fam,Site!='Duke' & Age!='endog'  & Type=='root'))))


save(RelAbund.withEndog.fam,file="rel_abund/RelAbund_withEndog_fam.RData")
####### Order level: Make dataframe of relative Order abundances in various subsets of samples -- with Endogenous plants #######
RelAbund.withEndog.ord<-data.frame("Taxon"=taxa_names(fullEco.ord),"Level"='ord') %>%
  mutate("RA_leaf5"=taxa_sums(subset_samples(fullEco.ord, Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord, Type=='leaf')))) %>%
  mutate("RA_root5"=taxa_sums(subset_samples(fullEco.ord, Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord, Type=='root')))) %>%
  mutate("RA_leaf3"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &  Type=='leaf')))) %>%
  mutate("RA_root3"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &  Type=='root')))) %>%
  mutate("RA_leaf_Jam"=taxa_sums(subset_samples(fullEco.ord,Site=='Jam' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Jam' &  Type=='leaf')))) %>%
  mutate("RA_root_Jam"=taxa_sums(subset_samples(fullEco.ord,Site=='Jam' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Jam' &  Type=='root')))) %>%
  mutate("RA_leaf_Mah"=taxa_sums(subset_samples(fullEco.ord,Site=='Mah' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Mah' &  Type=='leaf')))) %>%
  mutate("RA_root_Mah"=taxa_sums(subset_samples(fullEco.ord,Site=='Mah' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Mah' &  Type=='root')))) %>%
  mutate("RA_leaf_Sil"=taxa_sums(subset_samples(fullEco.ord,Site=='Sil' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Sil' &  Type=='leaf')))) %>%
  mutate("RA_root_Sil"=taxa_sums(subset_samples(fullEco.ord,Site=='Sil' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Sil' &  Type=='root')))) %>%
  mutate("RA_root_Par"=taxa_sums(subset_samples(fullEco.ord,Site=='Par' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Par' & Type=='root')))) %>%
  mutate("RA_root_Mil"=taxa_sums(subset_samples(fullEco.ord,Site=='Mil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Mil' & Type=='root')))) %>%
  mutate("RA_leaf_Par"=taxa_sums(subset_samples(fullEco.ord,Site=='Par' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Par' & Type=='leaf')))) %>%
  mutate("RA_leaf_Mil"=taxa_sums(subset_samples(fullEco.ord,Site=='Mil' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Mil' & Type=='leaf')))) %>%
  mutate("RA_root_GRHsoil"=taxa_sums(subset_samples(fullEco.ord,Treatment=='GRHsoil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Treatment=='GRHsoil' & Type=='root')))) %>%
  mutate("RA_soil_Jam"=taxa_sums(subset_samples(fullEco.ord,Site=='Jam' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Jam' & Type=='soil')))) %>%
  mutate("RA_soil_Mah"=taxa_sums(subset_samples(fullEco.ord,Site=='Mah' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Mah' & Type=='soil')))) %>%
  mutate("RA_soil_Sil"=taxa_sums(subset_samples(fullEco.ord,Site=='Sil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Sil' & Type=='soil')))) %>%
  mutate("RA_soil_Par"=taxa_sums(subset_samples(fullEco.ord,Site=='Par' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Par' & Type=='soil')))) %>%
  mutate("RA_soil_Mil"=taxa_sums(subset_samples(fullEco.ord,Site=='Mil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.ord,Site=='Mil' & Type=='soil')))) %>%
  mutate("RA_soil_GRH"=taxa_sums(subset_samples(fullEco.ord,Treatment=='GRHsoil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.ord,Treatment=='GRHsoil' & Type=='soil')))) %>%
  mutate("RA_leaf3_age2"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') & Age=='2' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age2"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') & Age=='2'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age3"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') & Age=='3' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age3"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') & Age=='3'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age4"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') & Age=='4' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age4"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') & Age=='4'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='root')))) %>%
  mutate("RA_leaf3_endog"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='leaf')))) %>%
  mutate("RA_root3_endog"=taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.ord,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='root'))))
save(RelAbund.withEndog.ord,file="rel_abund/RelAbund_withEndog_ord.RData")

####### Class level: Make dataframe of relative Class abundances in various subsets of samples -- with Endogenous plants #######
RelAbund.withEndog.cla<-data.frame("Taxon"=taxa_names(fullEco.cla),"Level"='cla') %>%
  mutate("RA_leaf5"=taxa_sums(subset_samples(fullEco.cla, Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla, Type=='leaf')))) %>%
  mutate("RA_root5"=taxa_sums(subset_samples(fullEco.cla, Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla, Type=='root')))) %>%
  mutate("RA_leaf3"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &  Type=='leaf')))) %>%
  mutate("RA_root3"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &  Type=='root')))) %>%
  mutate("RA_leaf_Jam"=taxa_sums(subset_samples(fullEco.cla,Site=='Jam' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Jam' &  Type=='leaf')))) %>%
  mutate("RA_root_Jam"=taxa_sums(subset_samples(fullEco.cla,Site=='Jam' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Jam' &  Type=='root')))) %>%
  mutate("RA_leaf_Mah"=taxa_sums(subset_samples(fullEco.cla,Site=='Mah' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Mah' &  Type=='leaf')))) %>%
  mutate("RA_root_Mah"=taxa_sums(subset_samples(fullEco.cla,Site=='Mah' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Mah' &  Type=='root')))) %>%
  mutate("RA_leaf_Sil"=taxa_sums(subset_samples(fullEco.cla,Site=='Sil' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Sil' &  Type=='leaf')))) %>%
  mutate("RA_root_Sil"=taxa_sums(subset_samples(fullEco.cla,Site=='Sil' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Sil' &  Type=='root')))) %>%
  mutate("RA_root_Par"=taxa_sums(subset_samples(fullEco.cla,Site=='Par' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Par' & Type=='root')))) %>%
  mutate("RA_root_Mil"=taxa_sums(subset_samples(fullEco.cla,Site=='Mil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Mil' & Type=='root')))) %>%
  mutate("RA_leaf_Par"=taxa_sums(subset_samples(fullEco.cla,Site=='Par' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Par' & Type=='leaf')))) %>%
  mutate("RA_leaf_Mil"=taxa_sums(subset_samples(fullEco.cla,Site=='Mil' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Mil' & Type=='leaf')))) %>%
  mutate("RA_root_GRHsoil"=taxa_sums(subset_samples(fullEco.cla,Treatment=='GRHsoil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Treatment=='GRHsoil' & Type=='root')))) %>%
  mutate("RA_soil_Jam"=taxa_sums(subset_samples(fullEco.cla,Site=='Jam' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Jam' & Type=='soil')))) %>%
  mutate("RA_soil_Mah"=taxa_sums(subset_samples(fullEco.cla,Site=='Mah' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Mah' & Type=='soil')))) %>%
  mutate("RA_soil_Sil"=taxa_sums(subset_samples(fullEco.cla,Site=='Sil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Sil' & Type=='soil')))) %>%
  mutate("RA_soil_Par"=taxa_sums(subset_samples(fullEco.cla,Site=='Par' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Par' & Type=='soil')))) %>%
  mutate("RA_soil_Mil"=taxa_sums(subset_samples(fullEco.cla,Site=='Mil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.cla,Site=='Mil' & Type=='soil')))) %>%
  mutate("RA_soil_GRH"=taxa_sums(subset_samples(fullEco.cla,Treatment=='GRHsoil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.cla,Treatment=='GRHsoil' & Type=='soil')))) %>%
  mutate("RA_leaf3_age2"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') & Age=='2' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age2"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') & Age=='2'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age3"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') & Age=='3' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age3"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') & Age=='3'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age4"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') & Age=='4' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age4"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') & Age=='4'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='root')))) %>%
  mutate("RA_leaf3_endog"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='leaf')))) %>%
  mutate("RA_root3_endog"=taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.cla,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='root'))))
save(RelAbund.withEndog.cla,file="rel_abund/RelAbund_withEndog_cla.RData")
####### Phylum level: Make dataframe of relative Phylum abundances in various subsets of samples -- with Endogenous plants #######
RelAbund.withEndog.phy<-data.frame("Taxon"=taxa_names(fullEco.phy),"Level"='phy') %>%
  mutate("RA_leaf5"=taxa_sums(subset_samples(fullEco.phy, Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy, Type=='leaf')))) %>%
  mutate("RA_root5"=taxa_sums(subset_samples(fullEco.phy, Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy, Type=='root')))) %>%
  mutate("RA_leaf3"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &  Type=='leaf')))) %>%
  mutate("RA_root3"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &  Type=='root')))) %>%
  mutate("RA_leaf_Jam"=taxa_sums(subset_samples(fullEco.phy,Site=='Jam' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Jam' &  Type=='leaf')))) %>%
  mutate("RA_root_Jam"=taxa_sums(subset_samples(fullEco.phy,Site=='Jam' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Jam' &  Type=='root')))) %>%
  mutate("RA_leaf_Mah"=taxa_sums(subset_samples(fullEco.phy,Site=='Mah' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Mah' &  Type=='leaf')))) %>%
  mutate("RA_root_Mah"=taxa_sums(subset_samples(fullEco.phy,Site=='Mah' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Mah' &  Type=='root')))) %>%
  mutate("RA_leaf_Sil"=taxa_sums(subset_samples(fullEco.phy,Site=='Sil' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Sil' &  Type=='leaf')))) %>%
  mutate("RA_root_Sil"=taxa_sums(subset_samples(fullEco.phy,Site=='Sil' &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Sil' &  Type=='root')))) %>%
  mutate("RA_root_Par"=taxa_sums(subset_samples(fullEco.phy,Site=='Par' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Par' & Type=='root')))) %>%
  mutate("RA_root_Mil"=taxa_sums(subset_samples(fullEco.phy,Site=='Mil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Mil' & Type=='root')))) %>%
  mutate("RA_leaf_Par"=taxa_sums(subset_samples(fullEco.phy,Site=='Par' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Par' & Type=='leaf')))) %>%
  mutate("RA_leaf_Mil"=taxa_sums(subset_samples(fullEco.phy,Site=='Mil' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Mil' & Type=='leaf')))) %>%
  mutate("RA_root_GRHsoil"=taxa_sums(subset_samples(fullEco.phy,Treatment=='GRHsoil' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Treatment=='GRHsoil' & Type=='root')))) %>%
  mutate("RA_soil_Jam"=taxa_sums(subset_samples(fullEco.phy,Site=='Jam' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Jam' & Type=='soil')))) %>%
  mutate("RA_soil_Mah"=taxa_sums(subset_samples(fullEco.phy,Site=='Mah' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Mah' & Type=='soil')))) %>%
  mutate("RA_soil_Sil"=taxa_sums(subset_samples(fullEco.phy,Site=='Sil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Sil' & Type=='soil')))) %>%
  mutate("RA_soil_Par"=taxa_sums(subset_samples(fullEco.phy,Site=='Par' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Par' & Type=='soil')))) %>%
  mutate("RA_soil_Mil"=taxa_sums(subset_samples(fullEco.phy,Site=='Mil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.phy,Site=='Mil' & Type=='soil')))) %>%
  mutate("RA_soil_GRH"=taxa_sums(subset_samples(fullEco.phy,Treatment=='GRHsoil' & Type=='soil'))/sum(taxa_sums(subset_samples(fullEco.phy,Treatment=='GRHsoil' & Type=='soil')))) %>%
  mutate("RA_leaf3_age2"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') & Age=='2' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age2"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') & Age=='2'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age3"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') & Age=='3' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age3"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') & Age=='3'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='root')))) %>%
  mutate("RA_leaf3_age4"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') & Age=='4' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='leaf')))) %>%
  mutate("RA_root3_age4"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') & Age=='4'  &  Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='root')))) %>%
  mutate("RA_leaf3_endog"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='leaf'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='leaf')))) %>%
  mutate("RA_root3_endog"=taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') & Age=='endog' & Type=='root'))/sum(taxa_sums(subset_samples(fullEco.phy,Site %in% c('Jam','Mah','Sil') &Age=='endog'  & Type=='root'))))
save(RelAbund.withEndog.phy,file="rel_abund/RelAbund_withEndog_phy.RData")

####### OTU level: Make dataframe of relative OTU abundances in various subsets of samples -- withOUT Endogenous plants #######
fieldEco.noEndog<-subset_samples(fieldEco,Genotype!='endog')
RelAbund.noEndog.otu<-data.frame("OTU_ID"=taxa_names(fieldEco.noEndog)) %>%
  mutate("allField"=taxa_sums(fieldEco.noEndog)/sum(taxa_sums(fieldEco.noEndog))) %>%
  mutate("leaf5"=taxa_sums(subset_samples(fieldEco.noEndog, Type=='leaf'))/sum(taxa_sums(subset_samples(fieldEco.noEndog, Type=='leaf')))) %>%
  mutate("root5"=taxa_sums(subset_samples(fieldEco.noEndog, Type=='root'))/sum(taxa_sums(subset_samples(fieldEco.noEndog, Type=='root')))) %>%
  mutate("leaf3"=taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &  Type=='leaf'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &  Type=='leaf')))) %>%
  mutate("root3"=taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &  Type=='root'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &  Type=='root')))) %>%
  mutate("leaf_Jam"=taxa_sums(subset_samples(fieldEco.noEndog,Site=='Jam' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site=='Jam' &  Type=='leaf')))) %>%
  mutate("root_Jam"=taxa_sums(subset_samples(fieldEco.noEndog,Site=='Jam' &  Type=='root'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site=='Jam' &  Type=='root')))) %>%
  mutate("leaf_Mah"=taxa_sums(subset_samples(fieldEco.noEndog,Site=='Mah' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site=='Mah' &  Type=='leaf')))) %>%
  mutate("root_Mah"=taxa_sums(subset_samples(fieldEco.noEndog,Site=='Mah' &  Type=='root'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site=='Mah' &  Type=='root')))) %>%
  mutate("leaf_Sil"=taxa_sums(subset_samples(fieldEco.noEndog,Site=='Sil' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site=='Sil' &  Type=='leaf')))) %>%
  mutate("root_Sil"=taxa_sums(subset_samples(fieldEco.noEndog,Site=='Sil' &  Type=='root'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site=='Sil' &  Type=='root')))) %>%
  mutate("leaf3_age2"=taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') & Age=='2' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='leaf')))) %>%
  mutate("root3_age2"=taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') & Age=='2'  &  Type=='root'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &Age=='2'  &  Type=='root')))) %>%
  mutate("leaf3_age3"=taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') & Age=='3' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='leaf')))) %>%
  mutate("root3_age3"=taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') & Age=='3'  &  Type=='root'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &Age=='3'  &  Type=='root')))) %>%
  mutate("leaf3_age4"=taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') & Age=='4' &  Type=='leaf'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='leaf')))) %>%
  mutate("root3_age4"=taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') & Age=='4'  &  Type=='root'))/sum(taxa_sums(subset_samples(fieldEco.noEndog,Site %in% c('Jam','Mah','Sil') &Age=='4'  &  Type=='root'))))
save(RelAbund.noEndog.otu,file="rel_abund/RelAbund_noEndog_otu.RData")

###### Get top 30 families at each site #######

top30fams.GRHvsField<-cbind(arrange(RelAbund.withEndog.fam,desc(RA_root5)) %>% mutate(RA_root5=RA_root5*100) %>% head(30) %>% select(Taxon,RA_root5) %>% plyr::rename(replace=c('Taxon'='all_field_roots')),
                            arrange(RelAbund.withEndog.fam,desc(RA_root_GRHsoil)) %>% mutate(RA_root_GRHsoil=RA_root_GRHsoil*100) %>% head(30) %>% select(Taxon,RA_root_GRHsoil) %>% plyr::rename(replace=c('Taxon'='GRHsoil_roots')),
                            arrange(RelAbund.withEndog.fam,desc(RA_soil_GRH)) %>% mutate(RA_soil_GRH=RA_soil_GRH*100) %>% head(30) %>% select(Taxon,RA_soil_GRH) %>% plyr::rename(replace=c('Taxon'='potting soil')))
                          
write.table(format(top30fams.GRHvsField,digits=1),file="tables/Table_S6_top30fams_GRHvsField.txt",sep='\t',col.names=TRUE,row.names=FALSE)


top30fams.root<-cbind(arrange(RelAbund.withEndog.fam,desc(RA_root_Jam)) %>% mutate(RA_root_Jam=RA_root_Jam*100) %>% head(30) %>% select(Taxon,RA_root_Jam) %>% plyr::rename(replace=c('Taxon'='Jam')),
                      arrange(RelAbund.withEndog.fam,desc(RA_root_Mah)) %>% mutate(RA_root_Mah=RA_root_Mah*100) %>% head(30) %>% select(Taxon,RA_root_Mah) %>% plyr::rename(replace=c('Taxon'='Mah')),
                      arrange(RelAbund.withEndog.fam,desc(RA_root_Sil)) %>% mutate(RA_root_Sil=RA_root_Sil*100) %>% head(30) %>% select(Taxon,RA_root_Sil) %>% plyr::rename(replace=c('Taxon'='Sil')),
                      arrange(RelAbund.withEndog.fam,desc(RA_root_Mil)) %>% mutate(RA_root_Mil=RA_root_Mil*100) %>% head(30) %>% select(Taxon,RA_root_Mil) %>% plyr::rename(replace=c('Taxon'='Mil')),
                      arrange(RelAbund.withEndog.fam,desc(RA_root_Par)) %>% mutate(RA_root_Par=RA_root_Par*100) %>% head(30) %>% select(Taxon,RA_root_Par) %>% plyr::rename(replace=c('Taxon'='Par')),
                      arrange(RelAbund.withEndog.fam,desc(RA_root_GRHsoil)) %>% mutate(RA_root_GRHsoil=RA_root_GRHsoil*100) %>% head(30) %>% select(Taxon,RA_root_GRHsoil) %>% plyr::rename(replace=c('Taxon'='GRHsoil')))

top30fams.leaf<-cbind(arrange(RelAbund.withEndog.fam,desc(RA_leaf_Jam)) %>% mutate(RA_leaf_Jam=RA_leaf_Jam*100) %>% head(30) %>% select(Taxon,RA_leaf_Jam) %>% plyr::rename(replace=c('Taxon'='Jam')),
                      arrange(RelAbund.withEndog.fam,desc(RA_leaf_Mah)) %>% mutate(RA_leaf_Mah=RA_leaf_Mah*100) %>% head(30) %>% select(Taxon,RA_leaf_Mah) %>% plyr::rename(replace=c('Taxon'='Mah')),
                      arrange(RelAbund.withEndog.fam,desc(RA_leaf_Sil)) %>% mutate(RA_leaf_Sil=RA_leaf_Sil*100) %>% head(30) %>% select(Taxon,RA_leaf_Sil) %>% plyr::rename(replace=c('Taxon'='Sil')),
                      arrange(RelAbund.withEndog.fam,desc(RA_leaf_Mil)) %>% mutate(RA_leaf_Mil=RA_leaf_Mil*100) %>% head(30) %>% select(Taxon,RA_leaf_Mil) %>% plyr::rename(replace=c('Taxon'='Mil')),
                      arrange(RelAbund.withEndog.fam,desc(RA_leaf_Par)) %>% mutate(RA_leaf_Par=RA_leaf_Par*100) %>% head(30) %>% select(Taxon,RA_leaf_Par) %>% plyr::rename(replace=c('Taxon'='Par')))

top30fams.soil<-cbind(arrange(RelAbund.withEndog.fam,desc(RA_soil_Jam)) %>% mutate(RA_soil_Jam=RA_soil_Jam*100) %>% head(30) %>% select(Taxon,RA_soil_Jam) %>% plyr::rename(replace=c('Taxon'='Jam')),
                      arrange(RelAbund.withEndog.fam,desc(RA_soil_Mah)) %>% mutate(RA_soil_Mah=RA_soil_Mah*100) %>% head(30) %>% select(Taxon,RA_soil_Mah) %>% plyr::rename(replace=c('Taxon'='Mah')),
                      arrange(RelAbund.withEndog.fam,desc(RA_soil_Sil)) %>% mutate(RA_soil_Sil=RA_soil_Sil*100) %>% head(30) %>% select(Taxon,RA_soil_Sil) %>% plyr::rename(replace=c('Taxon'='Sil')),
                      arrange(RelAbund.withEndog.fam,desc(RA_soil_Mil)) %>% mutate(RA_soil_Mil=RA_soil_Mil*100) %>% head(30) %>% select(Taxon,RA_soil_Mil) %>% plyr::rename(replace=c('Taxon'='Mil')),
                      arrange(RelAbund.withEndog.fam,desc(RA_soil_Par)) %>% mutate(RA_soil_Par=RA_soil_Par*100) %>% head(30) %>% select(Taxon,RA_soil_Par) %>% plyr::rename(replace=c('Taxon'='Par')),
                      arrange(RelAbund.withEndog.fam,desc(RA_soil_GRH)) %>% mutate(RA_soil_GRH=RA_soil_GRH*100) %>% head(30) %>% select(Taxon,RA_soil_GRH) %>% plyr::rename(replace=c('Taxon'='GRH')))

save(top30fams.leaf,file="rel_abund/top30fams_leaf.RData")
save(top30fams.root,file="rel_abund/top30fams_root.RData")
save(top30fams.soil,file="rel_abund/top30fams_soil.RData")

write.table(format(top30fams.leaf,digits=1),file="tables/Table_S2_top30fams_leaf.txt",sep='\t',col.names=TRUE,row.names=FALSE)
write.table(format(top30fams.soil,digits=1),file="tables/Table_S3_top30fams_soil.txt",sep='\t',col.names=TRUE,row.names=FALSE)
write.table(format(top30fams.root,digits=1),file="tables/Table_S4_top30fams_root.txt",sep='\t',col.names=TRUE,row.names=FALSE)

####### Save image #######
save.image(paste0("rel_abund/image_",date()))