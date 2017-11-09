setwd('/home/lucas/Main/Illumina/seed-microbes/')

##### Read source file #####
source('R/chapter1_source.R')

##### Preprocess ######
# Import sequence table generated from dada2, to be used as otu_table in phyloseq
seqtab <- readRDS('dada2/16S_final/16S_seqtab.Rds')
# Rename samples names to include "Sample" in front of each sample number
rownames(seqtab) <- paste("Sample",rownames(seqtab),sep="")

# Import taxa table, sample data, tree and create phyloseq object
taxtab <- readRDS('dada2/16S_final/16S_silva_taxtab_spp.Rds')
sd <- read.csv("data/mapping/field_mapping.csv", colClasses = c(rep('factor', 29), rep('numeric',31)))
rownames(sd) <- sd$SampleID
tree <- readRDS('dada2/16S_final/16S_tree.Rds')

##### Generate raw Phyloseq Object #####
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               tax_table(taxtab),
               sample_data(sd),
               phy_tree(tree))

# Remove putative bacterial contaminants (see contaminants.R file for details)
contaminants <- readRDS("dada2/16S_final/contaminant_otus.Rds")
ps <- prune_taxa(taxa_names(ps)[!taxa_names(ps) %in% contaminants],ps)
# Remove undefined and Chloroplast sequences
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps <- subset_taxa(ps, !Class %in% "Chloroplast")
ps <- subset_taxa(ps, Kingdom %in% "Bacteria")
# Remove negative controls
ps <-subset_samples(ps, Sample_type != "DNA_control")
# Remove low coverage samples with less than 500 sequences
ps <- prune_samples(names(which(sample_sums(ps) >= 500)),ps)

##### Correction for copy number estimation #####
# Copy number correction scripts go here - example from Maggie Wagner
## Load copy number estimates:
# CNCest<-as.data.frame(read.table('raw_data/Eco_Field_copynumest_forR.txt',sep='\t',header=TRUE))
# # do for fullEco dataset:
# fullOTU.CNC<-as.data.frame(as(otu_table(fullEco.nobadOTUs.highcoverage.thresholded),'matrix')) # copy OTU table from thresholded Phyloseq object
# fullOTU.CNC<-merge(CNCest,fullOTU.CNC,by.y='row.names',by.x='OTU_ID') # merge copy number estimates with OTU table
# rownames(fullOTU.CNC)<-fullOTU.CNC$OTU_ID # restore OTU_IDs as row names
# fullOTU.CNC<-fullOTU.CNC[,3:dim(fullOTU.CNC)[2]]/fullOTU.CNC[,2] # divide OTU counts by corresponding copy number estimate (Col 1 = OTU_ID, Col 2 = copy number estimates)
# fullOTU.CNC<-data.matrix(ceiling(fullOTU.CNC)) # round OTU counts up to nearest integer

##### Add alpha diversity metrics to sample data #####
ps <- add_alpha_diversity(ps)
ps <- add_alpha_rarefaction(ps)

####### Remove taxa below designated threshold  #######
threshold1<-kOverA(5,A=0) # set threshold values (require at least 5 samples with >0 reads)
ps <-filter_taxa(ps,threshold1,TRUE)

##### Rename OTU's/RSV's to represent taxonomic assignments #####
ps <- rename_otus(ps, otu_map=TRUE)

# Rename OTUs/RSVs
ps <- rename_otus(ps, otu_map = TRUE)

# Normalize count data using DESeq
ps.ds <- DESeq_varstab(ps, ~Tissue + Experiment)


ps.osu <- subset_samples(ps, Farm == "OSU" | Sample_type == "Seed_stock") %>% prune_taxa(taxa_sums(.)>0,.)
ps.osu.field <- subset_samples(ps.osu, Sample_type != "Seed_stock") %>% prune_taxa(taxa_sums(.)>0,.)
ps.osu.plant <- subset_samples(ps.osu.field, Tissue == "Seed" | Tissue == "Crown") %>% prune_taxa(taxa_sums(.)>0,.)
ps.osu.crown <- subset_samples(ps.osu.field, Tissue == "Crown") %>% prune_taxa(taxa_sums(.)>0,.)
ps.osu.seed <- subset_samples(ps.osu.field, Tissue == "Seed") %>% prune_taxa(taxa_sums(.)>0,.)

ps.2013 <- subset_samples(ps, Experiment == 1) %>% prune_taxa(taxa_sums(.)>0,.)
ps.2013.field <- subset_samples(ps.2013, Experiment == 1 & Sample_type != "Seed_stock" & Tissue != "Silk") %>% prune_taxa(taxa_sums(.)>0,.)
ps.2013.plant <- subset_samples(ps.2013.field, Tissue == "Seed" | Tissue == "Crown") %>% prune_taxa(taxa_sums(.)>0,.)
ps.2013.crown <- subset_samples(ps.2013.field, Tissue == "Crown") %>% prune_taxa(taxa_sums(.)>0,.)
ps.2013.seed <- subset_samples(ps.2013.field, Tissue == "Seed") %>% prune_taxa(taxa_sums(.)>0,.)

ps.2014 <- subset_samples(ps, Experiment == 1) %>% prune_taxa(taxa_sums(.)>0,.)
ps.2014.field <- subset_samples(ps.2014, Experiment == 1 & Sample_type != "Seed_stock") %>% prune_taxa(taxa_sums(.)>0,.)
ps.2014.plant <- subset_samples(ps.2014.field, Tissue == "Seed" | Tissue == "Crown") %>% prune_taxa(taxa_sums(.)>0,.)
ps.2014.crown <- subset_samples(ps.2014.field, Tissue == "Crown") %>% prune_taxa(taxa_sums(.)>0,.)
ps.2014.seed <- subset_samples(ps.2014.field, Tissue == "Seed") %>% prune_taxa(taxa_sums(.)>0,.)


saveRDS(ps,file="ps_objects/16S_field.ps.Rds")
saveRDS(ps.ds,file="ps_objects/16S_field.ps.vst.Rds")

saveRDS(ps.osu,file="ps_objects/16S_field.ps.osu.Rds")
saveRDS(ps.osu.field,file="ps_objects/16S_field.ps.osu.field.Rds")
saveRDS(ps.osu.plant,file="ps_objects/16S_field.ps.osu.plant.Rds")
saveRDS(ps.osu.crown,file="ps_objects/16S_field.ps.osu.crown.Rds")
saveRDS(ps.osu.seed,file="ps_objects/16S_field.ps.osu.seed.Rds")

saveRDS(ps.2013,file="ps_objects/16S_field.ps.2013.Rds")
saveRDS(ps.2013.field,file="ps_objects/16S_field.ps.2013.field.Rds")
saveRDS(ps.2013.plant,file="ps_objects/16S_field.ps.2013.plant.Rds")
saveRDS(ps.2013.crown,file="ps_objects/16S_field.ps.2013.crown.Rds")
saveRDS(ps.2013.seed,file="ps_objects/16S_field.ps.2013.seed.Rds")

saveRDS(ps.2014,file="ps_objects/16S_field.ps.2014.Rds")
saveRDS(ps.2014.field,file="ps_objects/16S_field.ps.2014.field.Rds")
saveRDS(ps.2014.plant,file="ps_objects/16S_field.ps.2014.plant.Rds")
saveRDS(ps.2014.crown,file="ps_objects/16S_field.ps.2014.crown.Rds")
saveRDS(ps.2014.seed,file="ps_objects/16S_field.ps.2014.seed.Rds")
