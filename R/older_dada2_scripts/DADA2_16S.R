library(dada2); packageVersion("dada2")

run1 <- "/home/lucas/Main/Illumina/seed-microbes/dada2/run1/16S"
run2 <- "/home/lucas/Main/Illumina/seed-microbes/dada2/run2/16S"

#Set path to whichever run I want to use in the following pipeline
path <- run2
setwd(path)


fns <- list.files(path)
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 2)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#plotQualityProfile(fnFs[[3]])
#plotQualityProfile(fnRs[[3]])

# Make directory and filenames for the filtered fastqs

filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)){
  dir.create(filt_path)
}
  
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter.
#I opted not to use trunclen parameter in order to save shorter ITS sequences that otherwise would have been truncated too short
#I tested different maxEE scores, and maxEE=12 filtered out ~25-50% of sequences

#For run 1, I set maxEE to c(8,8) and for run2 I set it to c(12,12) because many sequences were being filtered. This could be an issue

for(i in seq_along(fnFs)) {
fastqPairedFilter(c(fnFs[i], fnRs[i]), fout=c(filtFs[i], filtRs[i]),
                   truncQ=2, truncLen=c(250,200), maxN=0, maxEE=c(12,12), rm.phix=TRUE,
                   compress=TRUE, verbose=TRUE)
}

############################################################3
# NEW DADA2 tutorial suggestion. Automatically uses 1 million reads to determine error rate, randomize option selects random sequences
# FORWARD
# Learn error rates using standard method
errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
plotErrors(errF, nominalQ=TRUE)
saveRDS(errF,file='filtered/errF.RDS')
#
# # REVERSE
# # Learn error rates using standard method
errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)
plotErrors(errR, nominalQ=TRUE)
saveRDS(errR,file='filtered/errR.RDS')

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

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE, minOverlap = 20, maxMismatch=16)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

#####

st1 <- readRDS('/home/lucas/Main/Illumina/seed-microbes/dada2/run1/16S/seqtab.rds')
st2 <- readRDS('/home/lucas/Main/Illumina/seed-microbes/dada2/run2/16S/seqtab.rds')


# Take a look at the sequence tables
dim(st1)
#run1 yields a seqtab of 1542 reads for 16S
# Remove chimeras
st1.nochimera <- removeBimeraDenovo(st1, method="consensus", multithread=TRUE)
dim(st1.nochimera)
# 891 sequences left from run1

dim(st2)
#run2 yields a seqtab of 4016 reads for 16S
# Remove chimeras
st2.nochimera <- removeBimeraDenovo(st2, method="consensus", multithread=TRUE)
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


