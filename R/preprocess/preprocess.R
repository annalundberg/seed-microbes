setwd("~/Main/Illumina/16S/combined_fastq/bowtie/combined_together/usearch_open_ref/blast_assigned_taxonomy/")

library(phyloseq)
library(ggplot2)
library("DESeq2")
library(vegan)
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
#import mapping file and merge mapping file into phyloseq object
sd <- import_qiime_sample_data("field_mapping.txt")
ps_16S= merge_phyloseq(ps_16S, sd)
ps_16S



ps_ITS = import_biom("ITS_otu_table.biom", refseqfilename = "ITS_rep_set.fasta", parseFunction=parse_taxonomy_greengenes)
tax_table(ps_ITS) <- tax_table(ps_ITS)[,2:8]
sample_data(ps_ITS) <- sd


#set phyloseq object to generic name for script

ps <- ps_16S


######################### Clean up OTUs and Samples used in downstream analysis, create relative abundance, presence/absence OTU tables ############




#remove OTUs in negative controls from OTUS in samples
#controls = otu_table(ps)[,c("Sample69","Sample70","Sample71","Sample72","Sample190","Sample191")]
#controls = filter_taxa(controls, function(x) sum(x) > 0.1, T)
#controlOTUS = row.names(controls)
#badcontrols = otu_table(ps)[,c("Sample192","Sample96","Sample235","Sample236","Sample269","Sample270","Sample271","Sample272")]

#remove samples with less than 100 sequences
ps <- prune_samples(names(which(sample_sums(ps) >= 100)),ps)
# remove negative controls (for now)
ps <-subset_samples(ps, Tissue != "DNA_control")
#Only keep taxa occuring more than three times
ps <- prune_taxa(taxa_sums(ps) > 3, ps)
# create a presence/absence OTU table
ps.pa <- transform_sample_counts(ps,function(x)1*(x>0))
# create relative abundance table
ps.ra <- transform_sample_counts(ps,function(x)x/sum(x))
#use p/a table to identify sample counts per OTU, filter out those not in more than one sample
ps <- prune_taxa(taxa_sums(ps.pa) > 1, ps)



########### Normalize count data with DESeq2#############################

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


#Convert to normalized counts
ps <- DESeq_varstab(ps)


#################### Combine OTU table to different taxonomic levels #################

ps.genus <- tax_glom(ps, taxrank = "Genus")
ps.family <- tax_glom(ps.genus, taxrank = "Family")
