# filter taxa in OTU table
ps_crowns <- subset_samples(ps, Tissue== "Crown")
ps_crowns.pa <- transform_sample_counts(ps_crowns,function(x)1*(x>0))



tax_table(ps_crowns) <- tax_table(filter_taxa(ps_crowns.pa,function(x) sum(x)>4,T))

ps_crowns <- prune_samples(names(which(sample_sums(ps_crowns) >= 100)),ps_crowns)

ps.ra <- transform_sample_counts(ps,function(x)x/sum(x))

otu.disease.cortest <- corr.test(t(otu_table(ps_crowns.ra)), 
                                 as.matrix(as.numeric((sample_data(ps_crowns.ra)$disease))), 
                                 use="pairwise", adjust="BY")


ps_crowns <- ps

# This function finds all taxa in a phyloseq object that correlate with a specified covariate
# Phyloseq objects can be counts, relative abundance, or presence/absence
# number of 
# p < 0.05
# 

otu_corr <- function(phyloseq, p = 0.05, otu_cutoff = 0.){
  otu_list <- data.frame(matrix(ncol=14))
  colnames(otu_list) <- c("OTU","OTU_avg_ra","Nsamples", "Coeff", "p_value", "Rsq", "Rsq_adj", "Kingdom","Phylum","Class", "Order","Family","Genus", "Species")
  oturow <- 1
  for(taxon in taxa_names(phyloseq)){
    otu_samples <- data.frame(otu_table(phyloseq)[taxon,])
    prop_otus <- length(which(otu_samples[,] > 0))/length(otu_samples)
    #if(prop_otus )
    if(sum(otu_samples)<1){  # in other words, if using relative abundance data
      otu_samples <- otu_samples[1,which(otu_samples[,] > 0)] # only use samples that have taxa present, to protect against sampling bias
    }
    model_table <- t(otu_samples)
    rownames(model_table)
    model_table <- merge(model_table, subset(sample_data(phyloseq), select = disease), by="row.names")
    model <- lm(model_table$disease ~ model_table[,2],model_table)
    pvalue <- coef(summary(model))[2,4] 
    if(pvalue < 0.05){
      #rownames(otu_list)[oturow] <- taxon
      row <- c(taxon,rowMeans(otu_samples),length(otu_samples),coef(summary(model))[2,1],pvalue,summary(model)[8],summary(model)[8],tax_table(phyloseq)[taxon])
      print(row)
      otu_list[oturow,] <- row
      oturow = oturow + 1
    }
  }
  return(otu_list)
}

##### Old working funciton
#otu_corr <- function(phyloseq){
#otu_list <- data.frame(matrix(ncol=14))
#colnames(otu_list) <- c("OTU","OTU_avg_ra","Nsamples", "Coeff", "p_value", "Rsq", "Rsq_adj", "Kingdom","Phylum","Class", "Order","Family","Genus", "Species")
#oturow <- 1
##for(taxon in taxa_names(phyloseq)){
#otu_samples <- data.frame(otu_table(phyloseq)[taxon,])
#otu_samples <- otu_samples[1,which(otu_samples[,] > 0)]
#model_table <- t(otu_samples)
#rownames(model_table)
#model_table <- merge(model_table, subset(sample_data(phyloseq), select = disease), by="row.names")
#model <- lm(model_table$disease ~ model_table[,2],model_table)
#pvalue <- coef(summary(model))[2,4]
#if(pvalue < 0.05){
#rownames(otu_list)[oturow] <- taxon
#row <- c(taxon,rowMeans(otu_samples),length(otu_samples),coef(summary(model))[2,1],pvalue,summary(model)[8],summary(model)[8],tax_table(phyloseq)[taxon])
#print(row)
#otu_list[oturow,] <- row
#oturow = oturow + 1
#}
#}
#return(otu_list)
#}




sd <- data.frame(sample_data(ps))
alpha_metrics <- colnames(sd[,60:71])

sd<-subset.data.frame(sd,Year==2014 & Farm=="OSU")
ps <- subset_samples(ps,Year==2014 & Farm=="OSU")

for(metric in alpha_metrics){
  print(summary(lm(sd[,metric] ~ sd$log_copies_Fus_all, data = sd)))
}
 




sum(otu_samples)

otu_corr_table <- otu_corr(ps_crowns.pa)

write.table(otu_corr_table, "disease_corr_crowns_16S_species.csv", sep = ",")






