setwd("~/Main/Illumina/16S/combined_fastq/pear100_500/usearch_picked_otus")

library(phyloseq)
library(ggplot2)
library("DESeq2")
library(vegan)
library(psych)
library(biom)

########## Explore alpha diversity metrics using boxplots ##############
p <- plot_richness(ps1, color = "Tissue", x = "Tissue", measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))
p <- p + geom_boxplot(aes(fill = Tissue), alpha=0.1)
plot(p)



add_alpha_diversity <- function(ps){
  alpha.diversity <- estimate_richness(ps)
  sd <- data.frame(sample_data(ps))
  sd_alpha <- merge(sd,alpha.diversity, by = 0)
  sd_alpha <- sd_alpha[,2:ncol(sd_alpha)] # drop Row.names variable that got put in there
  rownames(sd_alpha) <- sd_alpha$SampleID
  sample_data(ps) <- sd_alpha
  return(ps)
}



#############################

##Alpha diversity (example from http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#alpha_diversity)
add_alpha_rarefaction <- function(ps){
  min_lib <- min(sample_sums(ps))
  nsamp=nsamples(ps)
  trials = 100
  richness = matrix(nrow = nsamp, ncol= trials)
  row.names(richness) = sample_names(ps)
  evenness = matrix(nrow = nsamp, ncol= trials)
  row.names(evenness) = sample_names(ps)
  set.seed(3)
  for (i in 1:100) {
    # Subsample
    r <- rarefy_even_depth(ps, sample.size = min_lib, verbose = FALSE, replace = TRUE)
    # Calculate richness
    rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
    richness[ ,i] <- rich
    # Calculate evenness
    even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
    evenness[ ,i] <- even
  }
  
  # Create a new dataframe to hold the means and standard deviations of richness estimates
  mean <- apply(richness, 1, mean)
  sd <- apply(richness, 1, sd)
  rich_stats <- data.frame(mean, sd)
  colnames(rich_stats) <- c("Richness_raref","Richness_raref_sd")
  
  # Create a new dataframe to hold the means and standard deviations of evenness estimates
  mean <- apply(evenness, 1, mean)
  sd <- apply(evenness, 1, sd)
  even_stats <- data.frame(mean, sd)
  colnames(even_stats) <- c("Inverse_Simpson_raref","Inverse_Simpson_raref_sd")
  
  alpha.raref <- merge(rich_stats, even_stats, by=0)
  rownames(alpha.raref) <- alpha.raref$Row.names
  alpha.raref <- alpha.raref[,2:ncol(alpha.raref)]
  sd <- data.frame(sample_data(ps))
  sd_alpha <- merge(sd,alpha.raref,by=0)
  sd_alpha <- sd_alpha[,2:ncol(sd_alpha)]
  rownames(sd_alpha) <- sd_alpha$SampleID
  sample_data(ps) <- sd_alpha
  
  return(ps)
}


