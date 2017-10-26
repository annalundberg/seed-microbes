library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library("DESeq2")
library(vegan)
library(psych)
library(biom)

setwd('/home/lucas/Main/Illumina/seed-microbes/dada2/16S_final')

set.seed(100)

##### Preprocess ######
# Import sequence table generated from dada2, to be used as otu_table in phyloseq
seqtab <- readRDS('16S_seqtab.Rds')
# Rename samples names to include "Sample" in front of each sample number
rownames(seqtab) <- paste("Sample",rownames(seqtab),sep="")

# Import taxa table, sample data, tree and create phyloseq object
taxtab <- readRDS('16S_silva_taxtab_spp.Rds')
sd <- read.csv("../../data/mapping/field_mapping.csv", colClasses = c(rep('factor', 29), rep('numeric',31)))
rownames(sd) <- sd$SampleID
tree <- readRDS('16S_tree.Rds')

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               tax_table(taxtab),
               sample_data(sd),
               phy_tree(tree))

# Remove putative bacterial contaminants (see contaminants.R file for details)
contaminants <- readRDS("contaminant_otus.Rds")
ps <- prune_taxa(taxa_names(ps)[!taxa_names(ps) %in% contaminants],ps)
# Remove undefined and Chloroplast sequences
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps <- subset_taxa(ps, !Class %in% "Chloroplast")


# # Compute prevalence of each feature, store as data.frame
# # Sourced from https://f1000research.com/articles/5-1492/v2
# prevdf = apply(X = otu_table(ps),
#                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
#                FUN = function(x){sum(x > 0)})
# # Add taxonomy and total read counts to this data.frame
# prevdf = data.frame(Prevalence = prevdf,
#                     TotalAbundance = taxa_sums(ps),
#                     tax_table(ps))
# # Subset to the remaining phyla
# prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
# ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
#   # Include a guess for parameter
#   geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
#   scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
#   facet_wrap(~Phylum) + theme(legend.position="none")

# remove negative controls
ps <-subset_samples(ps, Sample_type != "DNA_control")
#Only keep taxa occuring more than three times total
ps<- prune_taxa(taxa_sums(ps) > 3, ps)
# Create presence/absernce ps object
ps.pa <- transform_sample_counts(ps,function(x)1*(x>0))
# #Use p/a table to identify sample counts per OTU, filter out those not in more than one sample
ps <- prune_taxa(taxa_sums(ps.pa) > 1, ps)
# remove samples with less than 100 sequences
ps <- prune_samples(names(which(sample_sums(ps) >= 500)),ps)

# Create relative abundance ps object
ps.ra <- transform_sample_counts(ps,function(x)x/sum(x))

# How does the abundance distribution look?
qplot(sample_sums(ps)) + xlab("Counts per Sample")
qplot(log10(sample_sums(ps))) + xlab("Log10 Counts per Sample")

# Normalize count data using DESeq
DESeq_varstab <- function(phyloseq, model) {
  # phyloseq = the input phylose object that you want to get DESeq transformed counts for
  # design_variable = the design for the conversion to the DESeq object. must be in the form "as a function of", for example "~Host_Genus", must be a variable in the phyloseq object
  # Set variables to NULL
  deseq.vst = NULL
  geo_Means = NULL
  phyloseq.DESeq = NULL
  # Convert to a DESeq object
  deseq = phyloseq_to_deseq2(phyloseq, model)
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geo_Means = apply(counts(deseq), 1, gm_mean)
  # Check to see if any columns (samples) don't have any OTUs in them:
  if(sum(colSums(counts(deseq)) == 0) == 0) { # if all samples have taxa, go on
    # Now we step through the size factors, dispersions, and varience stabilization:
    deseq = estimateSizeFactors(deseq, geoMeans = geo_Means)
    deseq = estimateDispersions(deseq) # long step
    #plotDispEsts(deseq)
    deseq.vst = getVarianceStabilizedData(deseq)
    # replace negatives with zeros
    deseq.vst[deseq.vst <0] <- 0
    # add the varience stabilized otu numbers into the dataset:
    otu_table(phyloseq) <- otu_table(deseq.vst, taxa_are_rows = TRUE)
    # create a new object for the varience stabalized set
    phyloseq -> phyloseq.DESeq
    # And, filter any taxa that became 0s all the way across
    phyloseq.DESeq = filter_taxa(phyloseq.DESeq, function(x) sum(x) > 0.1, T)
    # return the new phyloseq object
    return(phyloseq.DESeq)
  } # end of IF loop 
  else {return("Error: your phyloseq object has samples with no taxa present.")}
} # end function


ps.ds <- DESeq_varstab(ps, ~1)

# Create a ps object with a simple log10 transformation
ps.log <- transform_sample_counts(ps, function(x) log(1 + x))

# Create presence/absernce ps object
ps.pa <- transform_sample_counts(ps,function(x)1*(x>0))


saveRDS(ps.ds,file="16S.ps.ds.Rds")
saveRDS(ps.pa,file="16S.ps.pa.Rds")
saveRDS(ps.ra,file="16S.ps.ra.Rds")
saveRDS(ps,file="16S.ps.Rds")
saveRDS(ps.log,file="16S.ps.log.Rds")
