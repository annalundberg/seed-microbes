
# In this script, I try pooling different amounts of samples for sequence variant estimation.
# As a result, I have generated sequence tables from different pooling sizes
# and a spreadsheet of sequence variant accumulation curves

library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")

amplicon = "16S"
setwd("/home/lucas/Main/Illumina/seed-microbes/dada2/")

path0 <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/combined/",amplicon,"/", sep="")
path1 <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/run1/",amplicon,"/", sep="")
path2 <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/run2/",amplicon,"/", sep="")


# Indentify .fastq sample files to be filtered
fns0 <- list.files(path0)
fns1 <- list.files(path1)
fns2 <- list.files(path2)
fns0 <- fns0[grepl('.fastq$', fns0)]
fns1 <- fns1[grepl('.fastq$', fns1)]
fns2 <- fns2[grepl('.fastq$', fns2)]
# Extract sample names from files to be filtered
sample.names0 <- sapply(strsplit(fns0, "_"), `[`, 2)
sample.names1 <- sapply(strsplit(fns1, "_"), `[`, 2)
sample.names2 <- sapply(strsplit(fns2, "_"), `[`, 2)
#Specify full path of files to be filtered
fns0 <- paste0(path0, fns0)
fns1 <- paste0(path1, fns1)
fns2 <- paste0(path2, fns2)


# from Run2, sample 17
# plotQualityProfile(fns1[[5]])
# plotQualityProfile(fns2[[9]])
# 
# # From Run1, sample 17
# plotQualityProfile(fnFs1[[1]])
# plotQualityProfile(fnRs1[[1]])

# Run1
# Create folders for filtered files
filtpath0 <- paste(path0, 'filtered', sep="")
filtpath1 <- paste(path1, 'filtered', sep="")
filtpath2 <- paste(path2, 'filtered', sep="") 
dir.create(filtpath0)
dir.create(filtpath1)
dir.create(filtpath2)
# Names of filtered files
filts0 <- file.path(filtpath0, paste0(amplicon,"_",sample.names0, "_filt.fastq.gz"))
filts1 <- file.path(filtpath1, paste0(amplicon,"_",sample.names1, "_filt.fastq.gz"))
filts2 <- file.path(filtpath2, paste0(amplicon,"_",sample.names2, "_filt.fastq.gz"))
names(filts0) <- sample.names0
names(filts1) <- sample.names1
names(filts2) <- sample.names2


# # Filter the pre-merged files
# # Combined
# for(i in seq_along(fns)){
#   fastqFilter(fns[i], filts[i], 
#               minLen=380,maxLen=460,maxN=0, maxEE=2, 
#               rm.phix=TRUE, compress=TRUE, verbose=TRUE)
# }
# 
# # Run1
# for(i in seq_along(fns1)){
#   fastqFilter(fns1[i], filts1[i], 
#               minLen=380,maxLen=460,maxN=0, maxEE=2, 
#               rm.phix=TRUE, compress=TRUE, verbose=TRUE)
# }
# # Run2
# for(i in seq_along(fns2)){
#   fastqFilter(fns2[i], filts2[i], 
#               minLen=380,maxLen=460,maxN=0, maxEE=2, 
#               rm.phix=TRUE, compress=TRUE, verbose=TRUE)
# }
# 


############################################################3
#####  Error Estimation

set.seed(100)

# COMBINED
# # Learn error rates using standard method
# err0 <- learnErrors(filts0, multithread=TRUE, randomize=TRUE, nreads = 3e+06)
# saveRDS(err0,file=paste0(path0,'err0.RDS'))
# 
# # RUN 1
# # Learn error rates using standard method
# err1 <- learnErrors(filts1, multithread=TRUE, randomize=TRUE, nreads = 3e+06)
# saveRDS(err1,file=paste0(path1,'err1.RDS'))
# #
# # RUN 2
# # Learn error rates using standard method
# err2 <- learnErrors(filts2, multithread=TRUE, randomize=TRUE, nreads = 3e+06)
# saveRDS(err2,file=paste0(path2,'err2.RDS'))


err0 <- readRDS(paste(path0,'err0.RDS', sep=''))
err1 <- readRDS(paste(path1,'err1.RDS', sep=''))
err2 <- readRDS(paste(path2,'err2.RDS', sep=''))

# plotErrors(err0, nominalQ=TRUE)
# plotErrors(err1, nominalQ=TRUE)
# plotErrors(err2, nominalQ=TRUE)

#### Create dataset to identify 

maxseqs.list <- c(1, 1e4, 1e5, 1e6, 1e7)
header <- c("Sample_set","Max_seqs","Iteration", "Pool_number","Pool_sample_names","Pool_numsamples", "Sequence_count", "Pool_variants", "New_variants", "Cum_variants")
df <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(df) <- header

sample.sets <- c("combined","run1","run2")
filt.sets <- list(filts0, filts1, filts2)
err.sets <- list(err0,err1,err2)
path.sets <- c(path0,path1,path2)


for(count in 1:length(sample.sets)){
  sample.set <- sample.sets[count]
  err <- err.sets[[count]]
  filts <- filt.sets[[count]]
  path <- path.sets[count]
  for(maxseqs in maxseqs.list){
    iters <- 3
    for(iter in 1:iters){
      #randomize samples
      rand.sample.names <- sample(names(filts))
      dds <- vector("list", length(rand.sample.names))
      names(dds) <- rand.sample.names
      #set variables
      drps <- vector("list")
      seqcount = 0
      pool.numsamples = 0
      sam.index = 0
      prev.sam.index = 0
      pool.ID = 0
      # go through each sample, pooling based on maxseqs 
      for(sam in rand.sample.names) {
        sam.index <- sam.index + 1
        pool.numsamples = pool.numsamples + 1
        derep <- derepFastq(filts[[sam]])
        drps[[pool.numsamples]] <- derep
        names(drps)[pool.numsamples] <- sam
        seqcount <- length(derep$uniques) + seqcount
        if(seqcount >= maxseqs || sam.index == length(rand.sample.names)){
          pool.ID <- pool.ID + 1
          pool.indexes <- (prev.sam.index+1):sam.index
          pool.sample.names <- toString(rand.sample.names[pool.indexes],sep=", ")
          cat(paste0("Determining sequence variants in pool ",pool.ID,".\n"))
          dds <- dada(drps, err=err, multithread=TRUE, pool=TRUE)
          #remove pooled derep to free up memory
          rm(drps)
          rm(derep)
          drps <- vector("list", 0)
          seqtab <- makeSequenceTable(dds)
          pool.seqs <- ncol(seqtab)
          if(pool.ID == 1){
            prevseqs <- 0
            seqtabs <- seqtab
          }else{
            prevseqs <- ncol(seqtabs)
            seqtabs <- mergeSequenceTables(seqtabs, seqtab) 
          }
          cum.seqs <- ncol(seqtabs)
          pool.new.seqs <- cum.seqs - prevseqs
          cat(paste0(pool.new.seqs," new sequence variants identified in pool; ", cum.seqs," cumulative\n"))
          row <- c(sample.set,maxseqs,iter,pool.ID, pool.sample.names, pool.numsamples,seqcount,pool.seqs,pool.new.seqs,cum.seqs)
          df[nrow(df)+1,] <- row
          seqcount = 0
          pool.numsamples = 0
          prev.sam.index <- sam.index
        }
      }
      write.csv(seqtabs, paste0(path,sample.set,"_seqtab_pooled_",maxseqs,"_",iter,".csv", sep = ",", col.names=NA))
      write.table(df, file=paste0(path,sample.set,"_accurve_pooled_",maxseqs,"_",iter,".tsv", quote=FALSE, sep='\t', col.names = NA))
    }
  }
}




