library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")

amplicon = "16S"
setwd("/home/lucas/Main/Illumina/seed-microbes/dada2/")

path <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/combined/",amplicon,"/", sep="")
path1 <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/run1/",amplicon,"/", sep="")
path2 <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/run2/",amplicon,"/", sep="")


# Indentify .fastq sample files to be filtered
fns <- list.files(path)
fns1 <- list.files(path1)
fns2 <- list.files(path2)
fns <- fns[grepl('.fastq$', fns)]
fns1 <- fns1[grepl('.fastq$', fns1)] # CHANGE if different file extensions
fns2 <- fns2[grepl('.fastq$', fns2)]
# Extract sample names from files to be filtered
sample.names <- sapply(strsplit(fns, "_"), `[`, 2)
sample.names1 <- sapply(strsplit(fns1, "_"), `[`, 2)
sample.names2 <- sapply(strsplit(fns2, "_"), `[`, 2)
#Specify full path of files to be filtered
fns <- paste0(path, fns)
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
filtpath <- paste(path, 'filtered', sep="")
filtpath1 <- paste(path1, 'filtered', sep="") ## filtered files go into the filtered/ subdirectory
filtpath2 <- paste(path2, 'filtered', sep="") ## filtered files go into the filtered/ subdirectory
dir.create(filtpath)
dir.create(filtpath1)
dir.create(filtpath2)
# Names of filtered files
filts <- file.path(filtpath, paste0(amplicon,"_",sample.names, "_filt.fastq.gz"))
filts1 <- file.path(filtpath1, paste0(amplicon,"_",sample.names1, "_filt.fastq.gz"))
filts2 <- file.path(filtpath2, paste0(amplicon,"_",sample.names2, "_filt.fastq.gz"))
names(filts) <- sample.names
names(filts1) <- sample.names1
names(filts2) <- sample.names2


# Filter the pre-merged files
# Combined
for(i in seq_along(fns)){
  fastqFilter(fns[i], filts[i], 
              minLen=380,maxLen=460,maxN=0, maxEE=2, 
              rm.phix=TRUE, compress=TRUE, verbose=TRUE)
}

# Run1
for(i in seq_along(fns1)){
  fastqFilter(fns1[i], filts1[i], 
              minLen=380,maxLen=460,maxN=0, maxEE=2, 
              rm.phix=TRUE, compress=TRUE, verbose=TRUE)
}
# Run2
for(i in seq_along(fns2)){
  fastqFilter(fns2[i], filts2[i], 
              minLen=380,maxLen=460,maxN=0, maxEE=2, 
              rm.phix=TRUE, compress=TRUE, verbose=TRUE)
}



############################################################3
#### Learn Errors ############

set.seed(100)

# COMBINED
# Learn error rates using standard method
err <- learnErrors(filts, multithread=TRUE, randomize=TRUE)
saveRDS(err,file=paste0(path,'err.RDS'))

# RUN 1
# Learn error rates using standard method
err1 <- learnErrors(filts1, multithread=TRUE, randomize=TRUE)
saveRDS(err1,file=paste0(path1,'err1.RDS'))
#
# RUN 2
# Learn error rates using standard method
err2 <- learnErrors(filts2, multithread=TRUE, randomize=TRUE)
saveRDS(err2,file=paste0(path2,'err2.RDS'))


err <- readRDS(paste(path,'err.RDS', sep=''))
err1 <- readRDS(paste(path1,'err1.RDS', sep=''))
err2 <- readRDS(paste(path2,'err2.RDS', sep=''))

plotErrors(err, nominalQ=TRUE)
plotErrors(err1, nominalQ=TRUE)
plotErrors(err2, nominalQ=TRUE)


dd$


####### Infer sequence variants #############
#### This script takes one sample at a time and builds an accumulation curve of newly identified sequence variants #### 
  
accurve1 = "combined"
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
counter = 0
for(sam in sample.names) {
  counter = counter + 1
  cat("Processing sample:", sam,",",counter,' of', length(sample.names), "samples \n")
  derep <- derepFastq(filts[[sam]])
  dd <- dada(derep, err=err, multithread=TRUE, pool=TRUE)
  dds1[[sam]] <- dd
  seqtab <- makeSequenceTable(dd)
  rownames(seqtab) <- sam
  if (counter<2){
    seqtabs <- seqtab
    prevseqs = 0
  }else{
    prevseqs <- ncol(seqtabs)
    seqtabs <- mergeSequenceTables(seqtabs, seqtab)
  }
  seqs <- ncol(seqtabs)
  sampleseqs <- ncol(seqtab)
  newseqs <- seqs - prevseqs
  cat(paste0(newseqs," new sequence variants identified in sample (",format(round(100*newseqs/sampleseqs, 0), nsmall = 0),"% of total found); ", seqs," cumulative\n"))
  accurve <- append(accurve,newseqs)
}
saveRDS(dds1, paste0(path1,"dada_run1.Rds"))
saveRDS(seqtabs1, paste0(path1,"seqtab_run1.Rds"))




accurve1 = "run1"
dds1 <- vector("list", length(sample.names1))
names(dds1) <- sample.names1
counter = 0
for(sam in sample.names1) {
  counter = counter + 1
  cat("Processing sample:", sam,",",counter,' of', length(sample.names1), "samples \n")
  derep <- derepFastq(filts1[[sam]])
  dd <- dada(derep, err=err1, multithread=TRUE, pool=TRUE)
  dds1[[sam]] <- dd
  seqtab <- makeSequenceTable(dd)
  rownames(seqtab) <- sam
  if (counter<2){
    seqtabs1 <- seqtab
    prevseqs = 0
  }else{
    prevseqs <- ncol(seqtabs1)
    seqtabs1 <- mergeSequenceTables(seqtabs1, seqtab)
  }
  seqs <- ncol(seqtabs1)
  sampleseqs <- ncol(seqtab)
  newseqs <- seqs - prevseqs
  cat(paste0(newseqs," new sequence variants identified in sample (",format(round(100*newseqs/sampleseqs, 0), nsmall = 0),"% of total found); ", seqs," cumulative\n"))
  accurve1 <- append(accurve1,newseqs)
}
saveRDS(dds1, paste0(path1,"dada_run1.Rds"))
saveRDS(seqtabs1, paste0(path1,"seqtab_run1.Rds"))


accurve2 = "run2"
dds2 <- vector("list", length(sample.names2))
names(dds2) <- sample.names2
counter = 0
for(sam in sample.names2) {
  counter = counter + 1
  cat("Processing sample:", sam,",",counter,' of', length(sample.names2), "samples \n")
  derep <- derepFastq(filts2[[sam]])
  dd <- dada(derep, err=err2, multithread=TRUE, pool=TRUE)
  dds2[[sam]] <- dd
  seqtab <- makeSequenceTable(dd)
  rownames(seqtab) <- sam
  if (counter<2){
    seqtabs2 <- seqtab
    prevseqs = 0
  }else{
    prevseqs <- ncol(seqtabs2)
    seqtabs2 <- mergeSequenceTables(seqtabs2, seqtab)
  }
  seqs <- ncol(seqtabs2)
  sampleseqs <- ncol(seqtab)
  newseqs <- seqs - prevseqs
  cat(paste0(newseqs," new sequence variants identified in sample (",format(round(100*newseqs/sampleseqs, 0), nsmall = 0),"% of total found); ", seqs," cumulative\n"))
  accurve2 <- append(accurve2,newseqs)
}
saveRDS(dds2, paste0(path2,"dada_run2.Rds"))
saveRDS(seqtabs2, paste0(path2,"seqtab_run2.Rds"))

st <- readRDS(paste0(path,"seqtab_combined.Rds"))
st1 <- readRDS(paste0(path1,"seqtab_run1.Rds"))
st2 <- readRDS(paste0(path2,"seqtab_run2.Rds"))


