library(dada2); packageVersion("dada2")

run1 <- "/home/lucas/Main/Illumina/seed-microbes/dada2/run1"
run2 <- "/home/lucas/Main/Illumina/seed-microbes/dada2/run2"

path <- run1
setwd(path)


fns <- list.files(run1)
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# plotQualityProfile(fnFs[[1]])
# plotQualityProfile(fnRs[[1]])

# Make directory and filenames for the filtered fastqs

filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) 
  dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter. 
#I opted not to use trunclen parameter in order to save shorter ITS sequences that otherwise would have been truncated too short
#I tested different maxEE scores, and maxEE=12 filtered out ~25-50% of sequences

# for(i in seq_along(fnFs)) {
#  fastqPairedFilter(c(fnFs[i], fnRs[i]), fout=c(filtFs[i], filtRs[i]),
#                    truncQ=2, maxN=0, maxEE=c(15,15), rm.phix=TRUE,
#                    compress=TRUE, verbose=TRUE)
# }

############################################################3
# NEW tutorial suggestion. Automatically uses 1 million reads to determine error rate, randomize option selects random sequences
# FORWARD
# # Learn error rates using standard method
# errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
# plotErrors(errF, nominalQ=TRUE)
# save(errF,file='filtered/errF.Rdata')
# 
# # REVERSE
# # Learn error rates using standard method
# errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)
# plotErrors(errR, nominalQ=TRUE)
# save(errR,file='filtered/errR.Rdata')

errF <- readRDS("/home/lucas/Main/Illumina/seed-microbes/dada2/run1/errF.RDS")
errR <- readRDS("/home/lucas/Main/Illumina/seed-microbes/dada2/run1/errR.RDS")

# Sample inference - all together, without subsetting - too much RAM

#derepF <- derepFastq(filtFs)
#derepR <- derepFastq(filtRs)
#add/remove pool=TRUE option to see if pooling samples results in more 
#dadaFs <- dada(derepF, err=errF, multithread=TRUE, pool=TRUE)
#dadaRs <- dada(derepR, err=errR, multithread=TRUE, pool=TRUE)

# Sample inference by subsetting - A quick script I wrote to subset & pool 32 samples at a time. 

lo <- c(1,33,65,97,129,161,193,225)
hi <- c(32, 64, 96, 128, 160, 192, 224, 256) 
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for (i in 1:8){
  print(i)
  # Sample inference
  derepF <- derepFastq(filtFs[lo[i]:hi[i]])
  derepR <- derepFastq(filtRs[lo[i]:hi[i]])
  #add/remove pool=TRUE option to see if pooling samples results in more 
  ddF <- dada(derepF, err=errF, multithread=TRUE, pool=TRUE)
  ddR <- dada(derepR, err=errR, multithread=TRUE, pool=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE, minOverlap = 20, maxMismatch=15)
  rm(derepF)
  rm(derepR)
  mergers[[lo[i]:hi[i]]] <- merger
}
dim(mergers)
seqtab <- makeSequenceTable(mergers)
dim(seqtable)
saveRDS(seqtab, paste(path,"seqtab_pooled32.rds", sep="/"))


mergers <- mergePairs(dadaFs, derepF, dadaRs, derepR, verbose=TRUE, minOverlap = 20, maxMismatch=12)
save(mergers, file=paste("mergers.Rdata",sep='/')) # CHANGE ME to where you want sequence table saved
rm(derepF);rm(derepR)
gc()
# seqtab <- makeSequenceTable(mergers)
# dim(seqtab)
# table(nchar(getSequences(seqtab)))
# View(seqtab)

#### Recommended pipeline for big data. Added parameters to mergePairs function

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE, minOverlap = 20, maxMismatch=15)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, paste(path,"seqtab.rds", sep="/"))

# Merge multiple runs (if necessary)
st1 <- readRDS("/path/to/run1/output/seqtab.rds")
st2 <- readRDS("/path/to/run2/output/seqtab.rds")
st3 <- readRDS("/path/to/run3/output/seqtab.rds")
st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/path/to/silva_nr_v123_train_set.fa.gz", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "/path/to/study/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "/pathto/study/tax_final.rds") # CHANGE ME ...

