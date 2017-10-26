library(ggplot2)
library(reshape2)
library(phyloseq)
library(dplyr)

setwd('/home/lucas/Main/Illumina/seed-microbes/dada2/16S_final')

# Import sequence table generated from dada2, to be used as otu_table in phyloseq
st12 <- readRDS('16S_seqtab.Rds')
# Rename samples names to include "Sample" in front of each sample number
rownames(st12) <- paste("Sample",rownames(st12),sep="")

# Import taxa table, sample data, and create phyloseq object
taxtab <- readRDS('16S_silva_taxtab_spp.Rds')
sd <- read.csv("../../data/mapping/field_mapping_with_plate3.csv")
rownames(sd) <- sd$SampleID
ps <- phyloseq(otu_table(st12, taxa_are_rows=FALSE),
               tax_table(taxtab),
               sample_data(sd))


# examine negative controls
ps.neg <-subset_samples(ps, Sample_type == "DNA_control")
ps.neg <- prune_taxa(taxa_sums(ps.neg) > 1, ps.neg)
sample_sums(ps.neg)



ps.neg.pa <- transform_sample_counts(ps.neg,function(x)1*(x>0))
sample_sums(ps.neg.pa) 
# How many species in each negative control?
# Sample185 Sample186 Sample187 Sample188 Sample189 Sample190 Sample191 Sample192 Sample235 Sample236 Sample269 Sample270 Sample271 
#        28        83        77        29       113        12         7       127        29        35        23        23        29 
# Sample272  Sample69  Sample70  Sample71  Sample72  Sample96 
#        35        45         3         3         5       479 

# Remove samples 185-189, as they are expected to have seed microbe contamination.
ps.neg.pa <- prune_samples(sample_names(ps.neg.pa)[6:19],ps.neg.pa)

# Find "core" taxa of negative controls
plot(taxa_sums(ps.neg.pa))
ps.neg.core <- prune_taxa(taxa_sums(ps.neg.pa) > 0.4*length(sample_names(ps.neg.pa)), ps.neg.pa)
plot(taxa_sums(ps.neg.core))
# Bacteria in over 40% of negative controls
tax_table(ps.neg.core)[,6]

# Identify "core" microbiome, i.e. taxa that exist in a high percentage of samples
ps <- prune_samples(names(which(sample_sums(ps) >= 100)),ps)
ps.pa <- transform_sample_counts(ps,function(x)1*(x>0))
plot(taxa_sums(ps.pa))
# Bacteria in over 50% of samples
ps.samples.core <- prune_taxa(taxa_sums(ps.pa) > 0.5*length(sample_names(ps.pa)), ps)
ps.samples.core <- prune_samples(sample_sums(ps.samples.core)>0,ps.samples.core)
ps.samples.core.ra <- transform_sample_counts(ps.samples.core,function(x)x/sum(x))

contaminants <- taxa_names(ps.neg.core)

#contaminants <- c(taxa_names(ps.samples.core),taxa_names(ps.neg.core))[duplicated(c(taxa_names(ps.samples.core),taxa_names(ps.neg.core)))]
tax_table(prune_taxa(contaminants,ps))[,6:7]

# A brief look in NCBI blast suggests either of these taxa could be endophytes or contaminants

# Propionibacterium acnes is a common skin bacterium , so I am judging it as a contaminant - even though it is also an endophyte of grapes.
# Acidovorax/Comamonas is often found associated with plants


ps.cont <- prune_taxa(contaminants,ps.samples.core.ra)


m.cont <- psmelt(ps.cont)
m.cont.mean <- aggregate(data=m.cont,Abundance~OTU+Sample_type+DNA_plate_number+Genus, FUN="mean")
colnames(m.cont.mean)[5] <- "Relative_Abundance"


genus_colors <- c(  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "red")

ggplot(data=m.cont.mean, aes(x = Sample_type, y = Relative_Abundance, fill = Genus))+ 
  facet_grid(~DNA_plate_number) +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5)) +
  geom_col(aes(x = Sample_type, y = Relative_Abundance, fill = Genus))


m <- psmelt(ps.samples.core.ra)
m.mean <- aggregate(data=m,Abundance~OTU+Tissue+DNA_plate_number+Genus, FUN="mean")
colnames(m.mean)[5] <- "Relative_Abundance"
ggplot(data=m.mean, aes(x = Tissue, y = Relative_Abundance, fill = Genus))+ 
  facet_grid(~DNA_plate_number) +
  scale_fill_manual(values = genus_colors) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5)) +
  geom_col()

# Decide to remove just the Propionibacterium "contaminant"
contaminants <- "AGATACCCTGGTAGTCCACGCTGTAAACGGTGGGTACTAGGTGTGGGGTCCATTCCACGGGTTCCGTGCCGTAGCTAACGCTTTAAGTACCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGGAATTGACGGGGCCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGTAGAACCTTACCTGGGTTTGACATGGATCGGGAGTGCTCAGAGATGGGTGTGCCTCTTTTGGGGTCGGTTCACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCACTGTTGCCAGCACGTTATGGTGGGGACTCAGTGGAGACCGCCGGGGTCAACTCGGAGGAAGGTGGGGATGAC"
saveRDS(contaminants,file="contaminant_otus.Rds")
