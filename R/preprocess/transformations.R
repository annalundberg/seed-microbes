#OTU Table transformations


taxtable <- tax

# remove OTU's that do not occur in more than one sample
which(taxa)


#check sample sum distribution
sums = data.frame(sum=sample_sums(ps))
ggplot(sums, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 100) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())


########
#Roo's DESeq script. Modified to make design= ~1 so we can run it without having to specify a design
DESeq_varstab <- function(phyloseq) {
  # phyloseq = the input phylose object that you want to get DESeq transformed counts for
  # design_variable = the design for the conversion to the DESeq object. must be in the form "as a function of", for example "~Host_Genus", must be a variable in the phyloseq object
  # Set variables to NULL
  deseq.vst = NULL
  geo_Means = NULL
  phyloseq.DESeq = NULL
  # Convert to a DESeq object
  deseq = phyloseq_to_deseq2(phyloseq, ~1)
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
    plotDispEsts(deseq)
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
#######################


#transform to presence/absence OTU table
ps.pa <- transform_sample_counts(ps,function(x)1*(x>0))
#transform to relative abundance table
ps.ra <- transform_sample_counts(ps,function(x)x/sum(x))



######################################################################################################################################
########## Adapted from script from Maggie Wagner to remove OTUs whose counts are less than 10% of the mean counts ###################
# Extract "non-rare" OTUs (count totals >10% of mean)
sum(sort(taxa_sums(ps.seeds),decreasing=TRUE)[1])/sum(taxa_sums(ps.seeds)) # top taxon is 8.3%
median(sort(taxa_sums(ps.seeds),decreasing=TRUE)) # 2
mean(taxa_sums(ps.seeds)) # 29.11898
mean(sort(taxa_sums(leaf3),decreasing=TRUE)[-1]) # without OTU_3, mean= 1081.96

com.ps.seed <-rownames(subset(data.frame("Count"=taxa_sums(ps.seed)),Count>0.1*mean(taxa_sums(ps.seed)))) # store "common" OTU_IDs: contains 1026 leaf OTUs
sort(taxa_sums(prune_taxa(com.ps.seed,ps.seed))/sum(taxa_sums(ps.seed)), decreasing=TRUE) # relative abundance (in entire dataset) 
sum(taxa_sums(prune_taxa(com.ps.seed,ps.seed)))/sum(taxa_sums(ps.seed)) # common set represents 96.46% of all observations

#Common OTUs, in over 10% of mean observations
com.ps.seed <- prune_taxa(taxa_sums(ps.seed)>0.1*mean(taxa_sums(ps.seed)), ps.seed)

###########################################################################################################################


