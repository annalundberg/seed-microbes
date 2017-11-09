# Source file for maize endophyte analysis

##### Color Palettes #####

palette <- colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

##### Load Packages ######
library(tidyr)
library(phyloseq)
library(genefilter)
library(plyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(lme4)
library(DESeq2)
library(reshape2)
library(VennDiagram)
library(mvabund)
library(boral)
library(genefilter)
library(wesanderson)

##### Zero-tolerant function to calculate grometric means #####
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

##### Add alpha-diversity metrics to sample data #####
add_alpha_diversity <- function(ps){
  alpha.diversity <- estimate_richness(ps)
  sd <- data.frame(sample_data(ps))
  sd_alpha <- merge(sd,alpha.diversity, by = 0)
  sd_alpha <- sd_alpha[,2:ncol(sd_alpha)] # drop Row.names variable that got put in there
  rownames(sd_alpha) <- sd_alpha$SampleID
  sample_data(ps) <- sd_alpha
  return(ps)
}

##### Add rarefaction alpha diversity metrics to sample data #####
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

##### Function to rename OTUs/RSVs using taxonomic assignment details #####
rename_otus <- function(ps, max.sp=2, otu_map = FALSE){
  taxtab <- tax_table(ps)
  taxa_list <- vector()
  seq_list <- vector()
  refseq <- character(length=0)
  for(row in 1:nrow(taxtab)){
    taxon <-taxtab[row,]
    seq_list <- c(seq_list,row.names(taxon))
    NAs <- sum(is.na(taxon))
    # Genus_species nomenclature for OTUs IDed to species
    if(NAs == 0){
      name <- paste(taxon[,6],taxon[,7], sep="_")
      # Restrict ambiguous species calls
      if(length(unlist(strsplit(name,"/"))) > max.sp){
        NAs <- 1
      }
    }
    # If otu is named to Genus, use "[Genus]_sp" convention
    if(NAs == 1){
      name <- paste(taxon[,6],"sp", sep="_")
    }
    # If OTU identified to higher taxonomic level, just use name of that level
    if(NAs > 1){
      name <- paste(taxon[,length(taxon)-NAs])
    }
    taxa_list <- c(taxa_list,name)
  }
  # Rename duplicate taxa in taxa_list, by adding numbers to the end
  for(taxon in 1:length(taxa_list)){
    matches <- taxa_list[taxon] == taxa_list
    nmatches <- sum(matches)
    if(nmatches > 1){
      taxa_list[matches] <- paste(taxa_list[matches],1:(nmatches),sep = "_")  
    }
  }
  if(otu_map == TRUE){
    for(i in 1:length(taxa_list)){
      refseq <- paste0(refseq,taxa_list[i], "\t",seq_list[i],"\n")
    }
    write.table(refseq,file="otu_name_mapfile.text")
  }
  taxa_names(ps) <- taxa_list
  return(ps)
}

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
  geo_Means = apply(counts(deseq), 1, gm_mean)
  # Check to see if any columns (samples) don't have any OTUs in them:
  if(sum(colSums(counts(deseq)) == 0) == 0) { # if all samples have taxa, go on
    # Now we step through the size factors, dispersions, and varience stabilization:
    deseq = estimateSizeFactors(deseq, geoMeans = geo_Means)
    deseq = estimateDispersions(deseq) # long step
    # #plotDispEsts(deseq)
    deseq.vst = getVarianceStabilizedData(deseq)
    # replace negatives with zeros
    deseq.vst[deseq.vst <0] <- 0
    # add the varience stabilized otu numbers into the dataset:
    otu_table(phyloseq) <- otu_table(deseq.vst.mat, taxa_are_rows = TRUE)
    # create a new object for the varience stabalized set
    phyloseq -> phyloseq.DESeq
    # And, filter any taxa that became 0s all the way across
    phyloseq.DESeq = filter_taxa(phyloseq.DESeq, function(x) sum(x) > 0.1, T)
    # return the new phyloseq object
    return(phyloseq.DESeq)
  } # end of IF loop 
  else {return("Error: your phyloseq object has samples with no taxa present.")}
} # end function

