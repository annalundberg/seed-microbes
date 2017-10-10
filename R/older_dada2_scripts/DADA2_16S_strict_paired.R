library(dada2); packageVersion("dada2")

path1 <- "/home/lucas/Main/Illumina/seed-microbes/dada2/run1/16S"
path2 <- "/home/lucas/Main/Illumina/seed-microbes/dada2/run2/16S"

#Set path to whichever run I want to use in the following pipeline

setwd("/home/lucas/Main/Illumina/seed-microbes/dada2")

# Run1
filtpathF1 <- file.path(path1, 'filteredF1') ## filtered files go into the filtered/ subdirectory
filtpathR1 <- file.path(path1, 'filteredR1') ## filtered files go into the filtered/ subdirectory

dir.create(filtpathF1)
dir.create(filtpathR1)

fns1 <- list.files(path1)
fastqs1 <- fns1[grepl('.fastq$', fns1)] # CHANGE if different file extensions

# Run2

filtpathF2 <- file.path(path2, 'filteredF2') ## filtered files go into the filtered/ subdirectory
filtpathR2 <- file.path(path2, 'filteredR2') ## filtered files go into the filtered/ subdirectory

dir.create(filtpathF2)
dir.create(filtpathR2)
fns2 <- list.files(path2)
fastqs2 <- fns2[grepl('.fastq$', fns2)] # CHANGE if different file extensions


# Run1
fnFs1 <- file.path(path1, fastqs1[grepl("_R1", fastqs1)] )# Just the forward read files
fnRs1 <- file.path(path1, fastqs1[grepl("_R2", fastqs1)] )# Just the reverse read files

sample.names1 <- sapply(strsplit(fastqs1[grepl("_R1", fastqs1)], "_"), `[`, 2)
sample.names.F1 <- paste(sample.names,"F", sep="_")
sample.names.R1 <- paste(sample.names,"R", sep="_")

# Run2
fnFs2 <- file.path(path2, fastqs2[grepl("_R1", fastqs2)] )# Just the forward read files
fnRs2 <- file.path(path2, fastqs2[grepl("_R2", fastqs2)] )# Just the reverse read files

sample.names2 <- sapply(strsplit(fastqs2[grepl("_R1", fastqs2)], "_"), `[`, 2)

sample.names.F2 <- paste(sample.names.F2,"F", sep="_")
sample.names.R2 <- paste(sample.names.R2,"R", sep="_")


# from Run2, sample 17
plotQualityProfile(fnFs2[[9]])
plotQualityProfile(fnRs2[[9]])

# From Run1, sample 17
plotQualityProfile(fnFs1[[1]])
plotQualityProfile(fnRs1[[1]])



# Run1 
filtFs1 <- file.path(filtpathF1, paste0(sample.names.F1, "_filt.fastq.gz"))
filtRs1 <- file.path(filtpathR1, paste0(sample.names.R1, "_filt.fastq.gz"))
# Run2
filtFs2 <- file.path(filtpathF2, paste0(sample.names.F2, "_filt.fastq.gz"))
filtRs2 <- file.path(filtpathR2, paste0(sample.names.R2, "_filt.fastq.gz"))


## TRUNCATION: Strict
# Run1
for(i in seq_along(fnFs1)) {
  fastqPairedFilter(c(fnFs1[i], fnRs1[i]), c(filtFs1[i], filtRs1[i]),
                    trimLeft=10,
                    truncLen=c(190,170), 
                    maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}

# Run2
for(i in seq_along(fnFs2)) {
  fastqPairedFilter(c(fnFs2[i], fnRs2[i]), c(filtFs2[i], filtRs2[i]),
                    trimLeft=10,
                    truncLen=c(190,170), 
                    maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}



############################################################3
# NEW DADA2 tutorial suggestion. Automatically uses 1 million reads to determine error rate, randomize option selects random sequences
## Run 1
# FORWARD
# Learn error rates using standard method
errF1 <- learnErrors(filtFs1, multithread=TRUE, randomize=TRUE)
plotErrors(errF1, nominalQ=TRUE)
saveRDS(errF1,file='errF1.RDS')
#
# # REVERSE
# # Learn error rates using standard method
errR1 <- learnErrors(filtRs1, multithread=TRUE, randomize=TRUE)
plotErrors(errR1, nominalQ=TRUE)
saveRDS(errR1,file='errR1.RDS')

## Run 2
# Learn error rates using standard method
errF2 <- learnErrors(filtFs2, multithread=TRUE, randomize=TRUE)
plotErrors(errF2, nominalQ=TRUE)
saveRDS(errF2,file='errF2.RDS')
#
# # REVERSE
# # Learn error rates using standard method
errR2 <- learnErrors(filtRs2, multithread=TRUE, randomize=TRUE)
plotErrors(errR2, nominalQ=TRUE)
saveRDS(errR2,file='errR2.RDS')




#errF <- readRDS("/home/lucas/Main/Illumina/seed-microbes/dada2/run1/16S/filtered/errF.RDS")
#errR <- readRDS("/home/lucas/Main/Illumina/seed-microbes/dada2/run1/16S/filtered/errR.RDS")

# Sample inference - all together, without subsetting - too much RAM

# derepF <- derepFastq(filtFs)
# # save(derepF,file='filtered/derepF.Rdata')
# derepR <- derepFastq(filtRs)
# # save(derepR,file='filtered/derepR.Rdata')
# 
# #derepF = load('filtered/derepF.Rdata')
# #derepR = load('filtered/derepR.Rdata')
# 
# #add/remove pool=TRUE option to see if pooling samples results in more 
# dadaFs <- dada(derepF, err=errF, multithread=TRUE, pool=TRUE)
# saveRDS(dadaFs,file='filtered/dadaFs.RDS')
# dadaRs <- dada(derepR, err=errR, multithread=TRUE, pool=TRUE)
# 
# 
# mergers <- mergePairs(dadaFs, derepF, dadaRs, derepR, verbose=TRUE, minOverlap = 20, maxMismatch=12)
# save(mergers, file=paste("mergers.Rdata",sep='/')) # CHANGE ME to where you want sequence table saved
# rm(derepF);rm(derepR)
# gc()
# seqtab <- makeSequenceTable(mergers)
# dim(seqtab)
# table(nchar(getSequences(seqtab)))
# View(seqtab)

#### Recommended pipeline for big data. Added parameters to mergePairs function, and have decent number of mergers

mergers1 <- vector("list", length(sample.names1))
names(mergers1) <- sample.names1
for(sam in 1:length(sample.names1)) {
  cat("Processing:", sam, "\n")
  derepF1 <- derepFastq(filtFs1[[sam]])
  ddF1 <- dada(derepF1, err=errF1, multithread=TRUE)
  derepR1 <- derepFastq(filtRs1[[sam]])
  ddR1 <- dada(derepR1, err=errR1, multithread=TRUE)
  merger1 <- mergePairs(ddF1, derepF1, ddR1, derepR1, verbose=TRUE, justConcatenate = TRUE)
  mergers1[[sam]] <- merger1
}

rm(derepF); rm(derepR)
seqtab1 <- makeSequenceTable(mergers1)
dim(seqtab1)


mergers2 <- vector("list", length(sample.names2))
names(mergers2) <- sample.names2
for(sam in 1:length(sample.names2)) {
  cat("Processing:", sam, "\n")
  derepF2 <- derepFastq(filtFs2[[sam]])
  ddF2 <- dada(derepF2, err=errF2, multithread=TRUE)
  derepR2 <- derepFastq(filtRs2[[sam]])
  ddR2 <- dada(derepR2, err=errR2, multithread=TRUE)
  merger2 <- mergePairs(ddF2, derepF2, ddR2, derepR2, verbose=TRUE, justConcatenate = TRUE)
  mergers2[[sam]] <- merger2
}

seqtab2 <- makeSequenceTable(mergers2)
dim(seqtab2)

# Construct sequence table and remove chimeras


#####

st1 <- readRDS('/home/lucas/Main/Illumina/seed-microbes/dada2/run1/16S/seqtab.rds')
st2 <- readRDS('/home/lucas/Main/Illumina/seed-microbes/dada2/run2/16S/seqtab.rds')


# Take a look at the sequence tables
dim(st1)
#run1 yields a seqtab of 1542 reads for 16S
# Remove chimeras
st1.nochimera <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE)
dim(st1.nochimera)
# 891 sequences left from run1

dim(st2)
#run2 yields a seqtab of 4016 reads for 16S
# Remove chimeras
st2.nochimera <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE)
dim(st2.nochimera)
# 1423 sequences are left from run2

## Merge sequence tables
st.all <- mergeSequenceTables(st1, st2)

# Error in mergeSequenceTables(st1, st2) : 
#   Duplicated sample names detected in the rownames.

totcol <- c(colnames(st1),colnames(st2))
length(totcol)
#5558
uniques <- totcol[!duplicated(totcol)]
length(uniques)
#5510
#this tells me there are only 548 shared sequence varaints between runs, about 10% of total variants between the two runs

totcol <- c(colnames(st1.nochimera),colnames(st2.nochimera))
length(totcol)
#2314
uniques <- totcol[!duplicated(totcol)]
length(uniques)
#2180 unique sequences, which means only 134 sequence variants shared between runs


#This one works to merge and add the sequence tables... at least it worked before
library(reshape2)
st1$names <- rownames(st1)
st2$names <- rownames(st2)
DF <- rbind(melt(st1,id="names"),melt(st2,id="names"))
st.all <- dcast(data=DF,formula=names ~ variable,fun.aggregate=sum)
dim(st.all)


