setwd("~/Main/Illumina/16S/combined_fastq/bowtie/combined_together/usearch_open_ref/blast_assigned_taxonomy/")

library(phyloseq)
library(ggplot2)
library("DESeq2")
library(vegan)
library(psych)
library(biom)





#ps_16S= import_biom("16S_otu_table.biom",
#                    refseqfilename = "16S_rep_set.fasta",
#                    treefilename = "16S_otus.tre",
#                    parseFunction=parse_taxonomy_greengenes)

#import simple biom file, add sample metadata
ps_16S= import_biom("16S_otu_table_simple.biom",
                    treefilename = "16S_otus.tre",
                    parseFunction=parse_taxonomy_greengenes)
# I removed the reference sequences for now [refseqfilename = "16S_rep_set.fasta"]
#import mapping file and merge mapping file into phyloseq object. I defined the column classes here. This might need to change.
sd <- data.frame(read.table("field_mapping.csv", sep=",", header=T, colClasses=c(rep('factor', 30), rep('numeric',30))))
row.names(sd) <- sd$SampleID
sample_data(ps_16S) <- sd




#ps_16S_DESeq2= import_biom("16S_otu_table_simple_DESeq2_filtered.biom",
#                    treefilename = "16S_otus.tre",
#                    parseFunction=parse_taxonomy_greengenes)
# left out: refseqfilename = "16S_rep_set.fasta"
#ps_16S_DESeq2 = merge_phyloseq(ps_16S_DESeq2,tax_table(ps_16S), sd)


#ps_16S_CSS= import_biom("16S_otu_table_simple_CSS_filtered.biom",
#                           treefilename = "16S_otus.tre",
#                           parseFunction=parse_taxonomy_greengenes)
# left out: refseqfilename = "16S_rep_set.fasta"
#ps_16S_CSS = merge_phyloseq(ps_16S_CSS,tax_table(ps_16S), sd)


ps_ITS = import_biom("ITS_otu_table.biom", parseFunction=parse_taxonomy_greengenes)
tax_table(ps_ITS) <- tax_table(ps_ITS)[,2:8]
#refseqfilename = "ITS_rep_set.fasta",
#add mapping file
sd <- data.frame(read.table("field_mapping.csv", sep=",", header=T, colClasses=c(rep('factor', 30), rep('numeric',30))))
row.names(sd) <- sd$SampleID
sample_data(ps_ITS) <- sd
ps_ITS = subset_taxa(ps_ITS, Kingdom=="Fungi")


################### set phyloseq object to generic name for script #################


ps <- ps_16S
ps <- ps_ITS


######################### Clean up OTUs and Samples used in downstream analysis, create relative abundance, presence/absence OTU tables ############



#remove OTUs in negative controls from OTUS in samples
#controls = otu_table(ps)[,c("Sample69","Sample70","Sample71","Sample72","Sample190","Sample191")]
#controls = filter_taxa(controls, function(x) sum(x) > 0.1, T)
#controlOTUS = row.names(controls)
#badcontrols = otu_table(ps)[,c("Sample192","Sample96","Sample235","Sample236","Sample269","Sample270","Sample271","Sample272")]


sums = data.frame(sum=sample_sums(ps))
median(sums[,1])
ggplot(sums, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth 16S") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#remove samples with less than 100 sequences
ps_16S <- prune_samples(names(which(sample_sums(ps_16S) >= 100)),ps_16S)
# remove negative controls (for now)
ps_16S<-subset_samples(ps_16S, Tissue != "DNA_control")
#Only keep taxa occuring more than three times
ps_16S<- prune_taxa(taxa_sums(ps_16S) > 3, ps_16S)
# create a presence/absence OTU table
ps_16S.pa <- transform_sample_counts(ps_16S,function(x)1*(x>0))
# create relative abundance table
ps_16S.ra <- transform_sample_counts(ps_16S,function(x)x/sum(x))
#use p/a table to identify sample counts per OTU, filter out those not in more than one sample
ps_16S<- prune_taxa(taxa_sums(ps_16S.pa) > 1, ps_16S)

sums = data.frame(sum=sample_sums(ps_ITS))
ggplot(sums, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 10) +
  ggtitle("Distribution of sample sequencing depth ITS") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#remove samples with less than 100 sequences
ps_ITS <- prune_samples(names(which(sample_sums(ps_ITS) >= 20)),ps_ITS)
# remove negative controls (for now)
ps_ITS <-subset_samples(ps_ITS, Tissue != "DNA_control")
#Only keep taxa occuring more than three times
ps_ITS <- prune_taxa(taxa_sums(ps_ITS) > 3, ps_ITS)
# create a presence/absence OTU table
ps_ITS.pa <- transform_sample_counts(ps_ITS,function(x)1*(x>0))
# create relative abundance table
ps_ITS.ra <- transform_sample_counts(ps_ITS,function(x)x/sum(x))
#use p/a table to identify sample counts per OTU, filter out those not in more than one sample
ps_ITS <- prune_taxa(taxa_sums(ps_ITS.pa) > 1, ps_ITS)

tax_table(ps_ITS)

########## Add alpha diversity data to sample_data, using functions defined in alpha_diversity.R ##############

ps <- add_alpha_diversity(ps)
ps <- add_alpha_rarefaction(ps)

ps <- ps_16S
write.table(sample_data(ps), "field_mapping_16Salphadiv.csv", sep=',',col.names=NA)

ps_ITS <- ps
########### Normalize count data with DESeq2#############################

#Convert to normalized counts, using a slighly modifed version of Roo's script, located in transformations.R
ps_16S <- DESeq_varstab(ps_16S)

ps_ITS <- DESeq_varstab(ps_ITS)


sort(sample_sums(ps))


#Do I want to rarefy? Rarefying to 1000 means I have to drop 24 samples that have low sequence abundance
ps<- rarefy_even_depth(ps,sample.size=1000,replace=F)

#################### Combine OTU table to different taxonomic levels #################

ps.genus <- tax_glom(ps, taxrank = "Genus")
ps.family <- tax_glom(ps.genus, taxrank = "Family")
ps.order <- tax_glom(ps.family, taxrank = "Order")


library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)

ps_16S_phylum$Phylum




ps_16S_class <- ps_16S %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class)

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)


ggplot(ps_16S_class, aes(x = Tissue, y = Abundance, fill = Class)) + 
  #facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = phylum_colors) +
  #scale_x_discrete(
  #  breaks = c("7/8", "8/4", "9/2", "10/6"),
   # labels = c("Jul", "Aug", "Sep", "Oct"), 
   # drop = FALSE
  #) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition of Corn and soil samples")


ps_ITS_class <- ps_ITS %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class)

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)


ggplot(ps_ITS_class, aes(x = Tissue, y = Abundance, fill = Class)) + 
  #facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = phylum_colors) +
  #scale_x_discrete(
  #  breaks = c("7/8", "8/4", "9/2", "10/6"),
  # labels = c("Jul", "Aug", "Sep", "Oct"), 
  # drop = FALSE
  #) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Class> 2%) \n") +
  ggtitle("Class Composition of Corn and soil samples")


########### Define all experiments and models to be used in ANOVA and PERMANOVA downstream data analysis ########################
########### Put all of them in a data_models list of experiments and corresponding models used ##################################

#sample data
sd <- data.frame(sample_data(ps))
sd

# Data from all years
F.full <- subset(sd, Experiment == 1)
F.full.model <- "response ~ Tissue+Year+Farm+Seed_source+Disinfection+Tissue:Year+Tissue:Farm+Tissue:Seed_source+Tissue:Disinfection+Year:Seed_source+Year:Disinfection+Farm:Seed_source+Farm:Disinfection+Seed_source:Disinfection+year_farm_section_block"
F.seeds <- subset(sd, Experiment ==1 & Tissue == "Seed")
F.seeds.model <- "response ~ Year+Farm+Seed_source+Disinfection+Year:Seed_source+Year:Disinfection+Farm:Seed_source+Farm:Disinfection+Seed_source:Disinfection+year_farm_section_block"
F.crowns <- subset(sd, Experiment ==1 & Tissue == "Crown")
F.crowns.model <- "response ~ Year+Farm+Seed_source+Disinfection+Year:Seed_source+Year:Disinfection+Farm:Seed_source+Farm:Disinfection+Seed_source:Disinfection+year_farm_section_block"

#Field data from 2013
F13.full <- subset(sd, Experiment ==1 & Year == 2013)
F13.full.model <- "response ~ Tissue+Seed_source+Disinfection+Tissue:Seed_source+Tissue:Disinfection+Seed_source:Disinfection+year_farm_section_block"
F13.seeds <- subset(sd, Experiment ==1 & Year == 2013 & Tissue == "Seed")
F13.seeds.model <- "response ~ Seed_source+Disinfection+Seed_source:Disinfection+year_farm_section_block"
F13.crowns <- subset(sd, Experiment ==1 & Year == 2013 & Tissue == "Crown")
F13.crowns.model <- "response ~ Seed_source+Disinfection+Seed_source:Disinfection+year_farm_section_block"

#OSU 2013,2014
F13.14.osu.full <- subset(sd, Experiment ==1 & Farm=="OSU")
F13.14.osu.full.model <- "response ~ Tissue+Year+Seed_source+Disinfection+Tissue:Year+Tissue:Seed_source+Tissue:Disinfection+Year:Seed_source+Year:Disinfection+Seed_source:Disinfection+year_farm_section_block"
F13.14.osu.seeds <- subset(sd, Experiment ==1 & Farm=="OSU" & Tissue == "Seed")
F13.14.osu.seeds.model <- "response ~ Year+Seed_source+Disinfection+Year:Seed_source+Year:Disinfection+Seed_source:Disinfection+year_farm_section_block"
F13.14.osu.crowns <- subset(sd, Experiment ==1 & Farm=="OSU" & Tissue == "Crown")
F13.14.osu.crowns.model <- "response ~ Year+Seed_source+Disinfection+Year:Seed_source+Year:Disinfection+Seed_source:Disinfection+year_farm_section_block"

#Field data from 2014
F14.full <- subset(sd, Experiment ==1 & Year == 2014)
F14.full.model <- "response ~ Tissue+Farm+Disinfection+Inoculation+Tissue:Farm+Tissue:Disinfection+Tissue:Inoculation+Farm:Disinfection+Farm:Inoculation+Disinfection:Inoculation+year_farm_section_block"
F14.seeds <- subset(sd, Experiment ==1 & Year == 2014 & Tissue == "Seed")
F14.seeds.model <- "response ~ Farm+Disinfection+Inoculation+Farm:Disinfection+Farm:Inoculation+Disinfection:Inoculation+year_farm_section_block"
F14.crowns <- subset(sd, Experiment ==1 & Year == 2014 & Tissue == "Crown")
F14.crowns.model <- "response ~ Farm+Disinfection+Inoculation+Farm:Disinfection+Farm:Inoculation+Disinfection:Inoculation+year_farm_section_block"

#OSU 2014
F14.osu.full <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "OSU")
F14.osu.full.model <- "response ~ Tissue+Seed_source+Disinfection+Inoculation+Tissue:Seed_source+Tissue:Disinfection+Tissue:Inoculation+Seed_source:Disinfection+Seed_source:Inoculation+Disinfection:Inoculation+year_farm_section_block"
F14.osu.seeds <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "OSU" & Tissue =="Seed")
F14.osu.seeds.model <- "response ~ Seed_source+Disinfection+Inoculation+Seed_source:Disinfection+Seed_source:Inoculation+Disinfection:Inoculation+year_farm_section_block"
F14.osu.crowns <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "OSU" & Tissue =="Crown")
F14.osu.crowns.model <- "response ~ Seed_source+Disinfection+Inoculation+Seed_source:Disinfection+Seed_source:Inoculation+Disinfection:Inoculation+year_farm_section_block"

# Adaptive Seeds
F14.as.full <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "AS")
F14.as.full.model <- "response ~ Tissue+Disinfection+Inoculation+Tissue:Disinfection+Tissue:Inoculation+Disinfection:Inoculation+year_farm_section_block"
F14.as.full.model <- "response ~ Tissue+Disinfection+Inoculation+Tissue:Disinfection+Tissue:Inoculation+Disinfection:Inoculation+year_farm_section_block"

F14.as.seeds <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "AS" & Tissue =="Seed")
F14.as.seeds.model <- "response ~ Disinfection+Inoculation+Disinfection:Inoculation+year_farm_section_block"
F14.as.crowns <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "AS" & Tissue =="Crown")
F14.as.crowns.model <- "response ~ Disinfection+Inoculation+Disinfection:Inoculation+year_farm_section_block"

#Pitchfork and Crow
F14.pc.full <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "PC")
F14.pc.full.model <- "response ~ Tissue+Disinfection+Inoculation+Tissue:Disinfection+Tissue:Inoculation+Disinfection:Inoculation+year_farm_section_block"
F14.pc.seeds <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "PC" & Tissue =="Seed")
F14.pc.seeds.model <- "response ~ Disinfection+Inoculation+Disinfection:Inoculation+year_farm_section_block"
F14.pc.crowns <- subset(sd, Experiment ==1 & Year == 2014 & Farm == "PC" & Tissue =="Crown")
F14.pc.crowns.model <-"response ~ Disinfection+Inoculation+Disinfection:Inoculation+year_farm_section_block"

experiment.names <- list("F.full","F.seeds","F.crowns","F13.full","F13.seeds","F13.crowns","F13.14.osu.full","F13.14.osu.seeds","F13.14.osu.crowns","F14.full","F14.seeds","F14.crowns","F14.osu.full","F14.osu.seeds","F14.osu.crowns","F14.as.full","F14.as.seeds","F14.as.crowns","F14.pc.full","F14.pc.seeds","F14.pc.crowns")
experiment.dfs <- list(F.full,F.seeds,F.crowns,F13.full,F13.seeds,F13.crowns,F13.14.osu.full,F13.14.osu.seeds,F13.14.osu.crowns,F14.full,F14.seeds,F14.crowns,F14.osu.full,F14.osu.seeds,F14.osu.crowns,F14.as.full,F14.as.seeds,F14.as.crowns,F14.pc.full,F14.pc.seeds,F14.pc.crowns)
models <- list(F.full.model,F.seeds.model,F.crowns.model,F13.full.model,F13.seeds.model,F13.crowns.model,F13.14.osu.full.model,F13.14.osu.seeds.model,F13.14.osu.crowns.model,F14.full.model,F14.seeds.model,F14.crowns.model,F14.osu.full.model,F14.osu.seeds.model,F14.osu.crowns.model,F14.as.full.model,F14.as.seeds.model,F14.as.crowns.model,F14.pc.full.model,F14.pc.seeds.model,F14.pc.crowns.model)
data_models <- list("experiment.names" = experiment.names, "experiment.dfs" = experiment.dfs, "models" = models)

mixed_data_models

######################################################################
##################################################################

#################################A list of distance matrices used #############################


#Distance matrices
bray_16S <- phyloseq::distance(ps_16S, method = "bray")
unifrac <- phyloseq::distance(ps, method = "unifrac")
wunifrac <- phyloseq::distance(ps, method = "wunifrac")
jaccard <- phyloseq::distance(ps, method ="jaccard")

bray.genus <- phyloseq::distance(ps.genus, method = "bray")
unifrac.genus <- phyloseq::distance(ps.genus, method = "unifrac")

bray.family <- phyloseq::distance(ps.family, method = "bray")
unifrac.family <- phyloseq::distance(ps.family, method = "unifrac")

################# Functions for Field Experiments ########################

# This function iterates over all experiments and models in a list of experiment/model combinations
# With the distance metric provided, it performs a permanova and outputs two tables, a p-table and a sum of squares table 
experiment_tables <- function(data_models, covariate=NULL, distance=NULL){
  ptable <- read.table(text = "",col.names = c("Tissue","Year","Farm","Seed_source","Disinfection","Inoculation","year_farm_section_block","Tissue:Year","Tissue:Farm","Tissue:Seed_source","Tissue:Disinfection","Tissue:Inoculation","Year:Seed_source","Year:Disinfection","Farm:Disinfection","Farm:Inoculation","Seed_source:Disinfection","Seed_source:Inoculation","Disinfection:Inoculation"))
  sstable <- read.table(text = "",col.names = c("Tissue","Year","Farm","Seed_source","Disinfection","Inoculation","year_farm_section_block","Tissue:Year","Tissue:Farm","Tissue:Seed_source","Tissue:Disinfection","Tissue:Inoculation","Year:Seed_source","Year:Disinfection","Farm:Disinfection","Farm:Inoculation","Seed_source:Disinfection","Seed_source:Inoculation","Disinfection:Inoculation","Residuals","Total","Percent_Explained"))
  for(i in 1:length(data_models$experiment.names)){
    print(paste("Running experiment ",i))
    experiment.name <- data_models$experiment.names[[i]]
    experiment.df <- data_models$experiment.dfs[[i]]
    model <- as.formula(data_models$models[[i]])
    # If the response is defined (dependent variable), run a multivariate anova
    if(!is.null(covariate)){
      response <- experiment.df[,covariate]
      result <- anova(lm(model,data=experiment.df))
      factors <- make.names(rownames(result))
      pvalues <- t(data.frame(result$`Pr(>F)`[1:(nrow(result)-1)]))
      colnames(pvalues) <- factors[1:(nrow(result)-1)]
      ssvalues <- t(data.frame(result$`Sum Sq`))
      colnames(ssvalues) <- factors
    }
    # If the distance matrix is defined, run an adonis
    if(!is.null(distance)){
      response <- as.dist(as.matrix(distance)[rownames(experiment.df),rownames(experiment.df)])
      result <- adonis(model,data=experiment.df, permutations=999)
      factors <- make.names(rownames(result$aov.tab[,]))
      pvalues <- t(data.frame(result$aov.tab$"Pr(>F)"[1:(length(factors)-2)]))
      colnames(pvalues) <- factors[1:(length(factors)-2)]
      ssvalues <- t(data.frame(result$aov.tab$"SumsOfSqs"[1:(length(factors))]))
      colnames(ssvalues) <- factors
    }
    # Populate each row the the table with the factors
    #ptable <- merge.data.frame(ptable,pvalues, all.y= T)[, union(names(ptable), names(pvalues))]
    for(col in 1:length(pvalues)){
      factor <- colnames(pvalues)[col]
      p <- pvalues[col]
      for(colname in colnames(ptable)){
        if(factor == colname){
          ptable[i,colname] <- p}}}
    rownames(ptable)[i] <- experiment.name
    #sstable <- merge.data.frame(sstable,ssvalues, all.y= T,sort=T)[, union(names(sstable), names(ssvalues))]
    for(col in 1:length(ssvalues)){
      factor <- colnames(ssvalues)[col]
      ss <- ssvalues[col]
      for(colname in colnames(sstable)){
        if(factor == colname){
          sstable[i,colname] <- ss}}}
    rownames(sstable)[i] <- experiment.name
    sstable$Percent_Explained[i] <- 100*(sstable$Total[i]-sstable$Residuals[i])/sstable$Total[i]
  }
  # Remove NA values
  ptable[is.na(ptable)] <- ""
  sstable[is.na(sstable)] <- ""
  return(list("ptable"=as.data.frame(ptable),"sstable"=as.data.frame(sstable)))
}




adonis_16S_bray <- experiment_tables(data_models, distance=bray)
adonis_16S_unifrac <- experiment_tables(data_models, distance=unifrac)
adonis_16S_jaccard <- experiment_tables(data_models, distance=jaccard)
adonis_16S_wunifrac <- experiment_tables(data_models, distance=wunifrac)

write.table(adonis_16S_bray$ptable,"stat_tables/16S_field_bray_ptable.csv", sep = ",", col.names=NA)
write.table(adonis_16S_bray$sstable,"stat_tables/16S_field_bray_sstable.csv", sep = ",", col.names=NA)

write.table(adonis_16S_unifrac$ptable,"stat_tables/16S_field_unifrac_ptable.csv", sep = ",", col.names=NA)
write.table(adonis_16S_unifrac$sstable,"stat_tables/16S_field_unifrac_sstable.csv", sep = ",", col.names=NA)

write.table(adonis_16S_jaccard$ptable,"stat_tables/16S_field_jaccard_ptable.csv", sep = ",", col.names=NA)
write.table(adonis_16S_jaccard$sstable,"stat_tables/16S_field_jaccard_sstable.csv", sep = ",", col.names=NA)

write.table(adonis_16S_wunifrac$ptable,"stat_tables/16S_field_wunifrac_ptable.csv", sep = ",", col.names=NA)
write.table(adonis_16S_wunifrac$sstable,"stat_tables/16S_field_wunifrac_sstable.csv", sep = ",", col.names=NA)

adonis_ITS_bray <- experiment_tables(data_models, distance=bray)
write.table(adonis_ITS_bray$ptable,"stat_tables/ITS_field_bray_ptable2.csv", sep = ",", col.names=NA)
write.table(adonis_ITS_bray$sstable,"stat_tables/ITS_field_bray_sstable.csv", sep = ",", col.names=NA)

adonis_ITS_jaccard <- experiment_tables(data_models, distance=jaccard)
write.table(adonis_ITS_jaccard$ptable,"stat_tables/ITS_field_jaccard_ptable.csv", sep = ",", col.names=NA)
write.table(adonis_ITS_jaccard$sstable,"stat_tables/ITS_field_jaccard_sstable.csv", sep = ",", col.names=NA)

field_simpsons_16S <- experiment_tables(data_models, covariate="Simpson")
write.table(field_simpsons_16S$ptable,"stat_tables/16S_field_simpson_ptable.csv", sep = ",", col.names=NA)
write.table(field_simpsons_16S$sstable,"stat_tables/16S_field_simpson_sstable.csv", sep = ",", col.names=NA)



################## Continuous Variables and Microbiome Structure ##################

experiment.names <- list("F.full","F.seeds","F.crowns","F13.full","F13.seeds","F13.crowns","F13.14.osu.full","F13.14.osu.seeds","F13.14.osu.crowns","F14.full","F14.seeds","F14.crowns","F14.osu.full","F14.osu.seeds","F14.osu.crowns","F14.as.full","F14.as.seeds","F14.as.crowns","F14.pc.full","F14.pc.seeds","F14.pc.crowns")
experiment.dfs <- list(F.full,F.seeds,F.crowns,F13.full,F13.seeds,F13.crowns,F13.14.osu.full,F13.14.osu.seeds,F13.14.osu.crowns,F14.full,F14.seeds,F14.crowns,F14.osu.full,F14.osu.seeds,F14.osu.crowns,F14.as.full,F14.as.seeds,F14.as.crowns,F14.pc.full,F14.pc.seeds,F14.pc.crowns)
models <- list(F.full.model,F.seeds.model,F.crowns.model,F13.full.model,F13.seeds.model,F13.crowns.model,F13.14.osu.full.model,F13.14.osu.seeds.model,F13.14.osu.crowns.model,F14.full.model,F14.seeds.model,F14.crowns.model,F14.osu.full.model,F14.osu.seeds.model,F14.osu.crowns.model,F14.as.full.model,F14.as.seeds.model,F14.as.crowns.model,F14.pc.full.model,F14.pc.seeds.model,F14.pc.crowns.model)
data_models <- list("experiment.names" = experiment.names, "experiment.dfs" = experiment.dfs, "models" = models)

mixed_data_models

experiment.df <- F14.full
distance <- bray
response <- as.dist(as.matrix(distance)[rownames(experiment.df),rownames(experiment.df)])
model <- as.formula("response ~ copies_Fus_all")
result <- adonis(model,data=experiment.df, permutations=999)



################### Microbe Sources ###################
F.sources <- subset(sd,Sample_type != "DNA_control", select=c(SampleID,DNA_sample_name,Experiment,Year,Farm,section,row,year_farm_section_block,Sample_type,Seed_source,Disinfection,Inoculation))
#F13.sources <- subset(sd, Year == 2013 & (Tissue == "Seed"|Tissue == "Crown"|Tissue == "Soil"))
#F14.sources <- subset(sd, Year==2014 & (Tissue == "Seed"|Tissue == "Crown"|Tissue == "Soil"))
#F13.14.osu.sources <- subset(sd, Experiment ==1 & Farm=="OSU" & (Tissue == "Seed"|Tissue == "Crown"|Tissue=="Soil"))
#F14.osu.sources <-subset(sd, Year==2014 & Farm =="OSU" & (Tissue == "Seed"|Tissue == "Crown"|Tissue == "Soil"))
#F14.as.sources <-subset(sd, Year==2014 & Farm =="AS" & (Tissue == "Seed"|Tissue == "Crown"|Tissue == "Soil"))  
#F14.pc.sources <-subset(sd, Year==2014 & Farm =="PC" & (Tissue == "Seed"|Tissue == "Crown"|Tissue == "Soil"))
View(F.sources)

Seed_src_2012as <- subset(F.sources,Sample_type == "Seed_stock" & Seed_source == 1)
Seed_src_2012cd <- subset(F.sources,Sample_type == "Seed_stock" & Seed_source == 2)
Seed_src_2013as <- subset(F.sources,Sample_type == "Seed_stock" & Seed_source == 3)
Seed_src_2013pc <- subset(F.sources,Sample_type == "Seed_stock" & Seed_source == 4)
Soil_source <- 
anova_
dist<- 

###############################3





# This function takes an experiment data frame, phyloseq object, and dependent variable
dist_model <- function(ps,sample_data,model,dep_var){
  dep <- sample_data$dep_var
  result <- aov(model,data=sample_data)
  return(result)
}

sample_data <- subset(F.full,is.na(F.full$disease)==F)
model <- F.full.model
dep_var <- "disease"
dep <- sample_data$disease



factors <- make.names(rownames(result$aov.tab[,]))
pvalues <- t(data.frame(result$aov.tab$"Pr(>F)"[1:(length(factors)-2)]))
colnames(pvalues) <- factors[1:(length(factors)-2)]

ssvalues <- t(data.frame(result$aov.tab$"SumsOfSqs"[1:(length(factors))]))
colnames(ssvalues) <- factors



# 
# And it computes
df_to_ps <- function(sd,ps){
  sdname <- deparse(substitute(sd))
  psname <- paste(sdname,".ps")
  F14.seeds.ps <- subset_samples(ps, SampleID %in% rownames(F14.seeds))
  #Only keep taxa occuring more than three times
  ps <- prune_ta)xa(taxa_sums(ps) > 3, ps)
  # create a presence/absence OTU table
  ps.pa <- transform_sample_counts(ps,function(x)1*(x>0))
  #use p/a table to identify sample counts per OTU, filter out those not in more than one sample
  ps <- prune_taxa(taxa_sums(ps.pa) > 1, ps)
}


betadisper <- 
anova()
permutest.betadisper(betadisper(bray_16S,sample_data(ps_16S)$Inoculation),pairwise=T,parallel=8,permutations=999)  


betadisper()

#Takes a df, factor and phyloseq object, and converts to a DESeq object and tests log2fold change with respect to the variable


ps <- ps_16S
F14.osu.seeds.ps <- ps
sample_data(F14.osu.seeds.ps) <- F14.osu.seeds
F14.osu.seeds.ps <- prune_taxa(taxa_sums(F.full.ps) > 0, F.full.ps)

deseq = phyloseq_to_deseq2(F14.osu.seeds.ps, ~Inoculation)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geo_Means = apply(counts(deseq), 1, gm_mean)
deseq <- estimateSizeFactors(deseq, geoMeans = geo_Means)
#deseq = estimateDispersions(deseq)
ps_DESeq = DESeq(deseq, test="Wald", fitType="parametric")

res = results(ps_DESeq, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

ord = ordinate(OFRF2014seeds.osu, method = "NMDS", distance = OFRF2014seeds.osu_bray)
plot_ordination(OFRF2014seeds.osu, ord, color = "Inoculation") + ggtitle("PCoA: Bray-Curtis")

#baseMean log2FoldChange     lfcSE     stat       pvalue         padj  Kingdom         Phylum               Class           Order
#3     59.59128       2.644583 0.5919291 4.467737 7.905152e-06 3.832418e-02 Bacteria Proteobacteria Gammaproteobacteria Xanthomonadales
#9028 312.55382       3.712706 0.5765149 6.439913 1.195419e-10 1.159079e-06 Bacteria Proteobacteria Gammaproteobacteria Xanthomonadales
#Family            Genus Species
#3    Xanthomonadaceae Stenotrophomonas    <NA>
#  9028 Xanthomonadaceae Stenotrophomonas    <NA>

inoc = subset_samples(OFRF2014seeds, Inoculation==1)
rowSums(otu_table(inoc)["9028",])
# 13964 individuals across all inoculated samples
inoc_ctl = subset_samples(OFRF2014seeds, Inoculation==0)
rowSums(otu_table(inoc_ctl)["9028",])
# 175 individuals across all control samples





############### Inoculation Effects 14.osu.seeds ################
ps <- ps_ITS

OFRF2014seeds.osu <- subset_samples(ps.genus, Year == 2014 & Experiment==1 & Farm=="AS" & Tissue == "Seed")
prune_taxa(taxa_sums(OFRF2014seeds.osu) > 0, OFRF2014seeds.osu)
ps_DESeq <- phyloseq_to_deseq2(ps_16S, ~Inoculation)
#
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geo_Means = apply(counts(ps_DESeq), 1, gm_mean)
ps_DESeq <- estimateSizeFactors(ps_DESeq, geoMeans = geo_Means)
#
ps_DESeq = DESeq(ps_DESeq, test="Wald", fitType="parametric")
res = results(ps_DESeq, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)

write.table(sigtab,"stat_tables/16S_foldchange_inoc_14_osu_seed_family.csv", sep=',',col.names = NA)

inoc = subset_samples(OFRF2014seeds.osu, Inoculation==1)
rowSums(otu_table(inoc)["41",])
rowMeans(otu_table(inoc)["41",])
# 2264 individuals across all inoculated samples
inoc_otu = otu_table(inoc)
inoc_otu[inoc_otu > 0] <- 1
dim(inoc_otu)
rowSums(inoc_otu["41",])
# 11/16 samples have Rhodococcus

inoc_ctl = subset_samples(OFRF2014seeds.osu, Inoculation==0)
rowSums(otu_table(inoc_ctl)["41",])
rowMeans(otu_table(inoc_ctl)["41",])
# 129 individuals across all control samples
inoc_ctl_otu = otu_table(inoc_ctl)
inoc_ctl_otu[inoc_ctl_otu > 0] <- 1
dim(inoc_ctl_otu)
rowSums(inoc_ctl_otu["41",])
# 7/16 samples have Rhodococcus


#Sort table by Family
sigtab <-sigtab[with(sigtab, order(Family,-log2FoldChange)), ]
sigtab[11,11]

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



p <- plot_composition()




############### Inoculation Effects 14.osu.crowns ################
ps <- ps_16S

OFRF2014 <- subset_samples(ps.family, Year == 2014 & Experiment==1 & Tissue == "Seed")
prune_taxa(taxa_sums(OFRF2014) > 0, OFRF2014)
ps_DESeq <- phyloseq_to_deseq2(OFRF2014, ~Inoculation)
#
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geo_Means = apply(counts(ps_DESeq), 1, gm_mean)
ps_DESeq <- estimateSizeFactors(ps_DESeq, geoMeans = geo_Means)
#
ps_DESeq = DESeq(ps_DESeq, test="Wald", fitType="parametric")
res = results(ps_DESeq, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_ITS)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)

write.table(sigtab,"stat_tables/ITS_foldchange_inoc_14_osu_seed_family.csv", sep=',',col.names = NA)

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
