library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")


amplicon = "16S"
setwd("/home/lucas/Main/Illumina/seed-microbes/dada2/pooling")

path0 <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/combined/",amplicon,"/pooling/", sep="")
path1 <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/run1/",amplicon,"/pooling/", sep="")
path2 <- paste("/home/lucas/Main/Illumina/seed-microbes/dada2/run2/",amplicon,"/pooling/", sep="")


# st0 <- read.csv("combined_seqtab_pooled_1e+07_1.csv", row.names=1)
# st1 <- read.csv("run1_seqtab_pooled_1e+06_1.csv", row.names=1)
# st2 <- read.csv("run2_seqtab_pooled_1e+06_1.csv", row.names=1)

st0 <- readRDS(paste0(path0,"seqtab_pooled_1e+07_1.Rds"))
st1 <- readRDS(paste0(path1,"seqtab_pooled_1e+06_1.Rds"))
st2 <- readRDS(paste0(path2,"seqtab_pooled_1e+06_1.Rds"))

ncol(st0) # 1783
ncol(st1) # 1228
ncol(st2) # 1223

# Remove Chimeras
st0_nc <- removeBimeraDenovo(st0, method="pooled", multithread=TRUE)
st1_nc <- removeBimeraDenovo(st1, method="pooled", multithread=TRUE)
st2_nc <- removeBimeraDenovo(st2, method="pooled", multithread=TRUE)

ncol(st0_nc) # 1242
100-100*(ncol(st0_nc)/ncol(st0)) # 30.34 % of sequence variants are chimeric
100-100*(sum(colSums(st0_nc))/sum(colSums(st0))) # 13.3% of all sequences are chimeric
ncol(st1_nc) # 782
100-100*(ncol(st1_nc)/ncol(st1)) # 36.32 of sequence variants are chimeric
100-100*(sum(colSums(st1_nc))/sum(colSums(st1))) # 14.5% of all sequences are chimeric
ncol(st2_nc) # 870
100-100*(ncol(st2_nc)/ncol(st2)) # 28.86 of sequence variants are chimeric
100-100*(sum(colSums(st2_nc))/sum(colSums(st2))) # 9.54% of all sequences are chimeric


# st0_chimeras <- st0[,!(colnames(st0) %in% colnames(st0_nc))]
# st1_chimeras <- st1[,!(colnames(st1) %in% colnames(st1_nc))]
# st2_chimeras <- st2[,!(colnames(st2) %in% colnames(st2_nc))]
# 
# taxa0_chimeras <- assignTaxonomy(st0_chimeras, "../../reference_database/RDP/rdp_train_set_16.fa.gz", multithread=TRUE)
# taxa1_chimeras <- assignTaxonomy(st1_chimeras, "../../reference_database/RDP/rdp_train_set_16.fa.gz", multithread=TRUE)
# taxa2_chimeras <- assignTaxonomy(st2_chimeras, "../../reference_database/RDP/rdp_train_set_16.fa.gz", multithread=TRUE)
# 


# I feel good about these numbers - I'll go ahead and go with the chimera-free sequence tables

st0 <- st0_nc
st1 <- st1_nc
st2 <- st2_nc



# Now, to look into how to combine sequencing runs 1 and 2

# How similar are runs 1 and 2? Let's start by removing samples from run2 that are not in run 1.
st2r <- st2[rownames(st1),]
dim(st2r) # 256 870
# remove any taxa that now have zero sums
st2r <- st2r[,colSums(st2r)!=0]
dim(st2r) # 256 856 ... 14 zero-abundance taxa were removed from st2r

seqs1 <- colnames(st1)
seqs2r <- colnames(st2r)

seqs_1.2r <- c(seqs1,seqs2r)
length(seqs_1.2r) # 1647
length(unique(seqs_1.2r)) # 1151

# Which sequences are shared between run1 and run2?
shared_1.2r <- seqs_1.2r[duplicated(seqs_1.2r)]
length(shared_1.2r) # 496 are shared sequence variants between runs

# How many overall sequences does this represent?

sum(colSums(st1[,shared_1.2r])) # 359782 shared total
100*(sum(colSums(st1[,shared_1.2r]))/sum(colSums(st1))) # 87.1% of sequences in run1 are shared with run2

sum(colSums(st2[,shared_1.2r])) # 754870 shared total
100*(sum(colSums(st2[,shared_1.2r]))/sum(colSums(st2))) # 88.6% of sequences in run2 are shared with run1
# 84% of all sequences in run1 are shared with run2

# Ideally, there will be no highly abundant taxa in either run that are unique to that run

unique_st1 <- st1[,!(colnames(st1) %in% shared_1.2r)] # Top 10 sequence abundances: 2549, 2429, 1179, 1127, 973, 929
max(colSums(unique_st1)) # 2549
100*max(colSums(unique_st1))/sum(colSums(st1)) # Most abundant unique variant is 0.62% of all sequences in run1, 6 of 1000
unique_st2r <- st2r[,!(colnames(st2r) %in% shared_1.2r)] 
max(colSums(unique_st2r)) # 3521
100*max(colSums(unique_st2r))/sum(colSums(st2)) # Most abundant unique variant is 0.41% of all sequences in run2, 4 of 1000

# I feel fairly confident that I can merge sequence tables between run1 and run2... what if I merge before filtering and learn errors?

seqs0 <- colnames(st0)
seqs1 <- colnames(st1)
seqs2 <- colnames(st2)

seqs_1.2 <- unique(c(seqs1,seqs2))
length(seqs_1.2) # 1156
length(seqs0) # 1242

seqs_0_1.2 <- c(seqs_1.2,seqs0)
shared_0_1.2 <- seqs_0_1.2[duplicated(seqs_0_1.2)]
length(shared_0_1.2) # 1025 sequence variants are shared when combined
100*length(shared_0_1.2)/length(seqs0) # That's 82.5 % of sequences
100*sum(colSums(st0[,shared_0_1.2]))/sum(colSums(st0)) # 96.7 % of sequences are shared with combined run1, run2

# I would assume that the unshared taxa do not represent highly abundant taxa - rather, they would be rare variants that are picked up with more pooling

unique_st0 <- st0[,!(colnames(st0) %in% shared_0_1.2)]
max(colSums(unique_st0)) # 3826 sequences for the most common unique variant
100*max(colSums(unique_st0))/sum(colSums(st0)) # 0.31% of all sequences

#### Merge sequence tables from run1 and run2 - dada2 does not allow same samples to be in two separate sequence tables
library(reshape2)
# Convert sequence tables to data frames and add column representing row names
st1_df <- data.frame(st1)
st2_df <- data.frame(st2)
st1_df$names <- rownames(st1)
st2_df$names <- rownames(st2)
# Merge data frames, using sample names, and sum duplicate entries
DF <- rbind(melt(st1_df,id="names"),melt(st2_df,id="names"))
st12 <- dcast(data=DF,formula=names ~ variable,fun.aggregate=sum)
dim(st12) # 288 1157 
# remove "names" column from data frames
rownames(st12) <- st12$names
st12 <- subset(st12,select=-names)
# Convert back to a matrix
st12 <- data.matrix(st12,rownames.force=TRUE)
dim(st12) # 288 1156
sum(colSums(st12)) # 1264833
sum(colSums(st0)) # 1235719

# More sequences are retained when merging runs 1 and 2, even though there are fewer sequence variants than in the pre-combined run.
# I will go with st12 (i.e. combining runs 1 and 2 after error estimation and sequence variant inference)

##### Save finalized sequence table, and assign taxonomy ############

setwd(paste0("/home/lucas/Main/Illumina/seed-microbes/dada2/",amplicon,"_final"))
saveRDS(st12,file=paste0(amplicon,"_seqtab.Rds"))
write.csv(st12, paste0(amplicon,"_seqtab.csv", sep = ",", col.names=NA))


taxa12_silva <- assignTaxonomy(st12, "../../reference_database/SILVA/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa12_silva_spp <- addSpecies(taxa12_silva, "../../reference_database/SILVA/silva_species_assignment_v123.fa.gz", verbose=TRUE, allowMultiple=TRUE)
saveRDS(taxa12_silva,file=paste0(amplicon,"_silva_taxtab.Rds"))
write.csv(taxa12_silva, paste0(amplicon,"_silva_taxtab.csv", col.names=NA))
saveRDS(taxa12_silva_spp,file=paste0(amplicon,"_silva_taxtab_spp.Rds"))
write.csv(taxa12_silva_spp, paste0(amplicon,"_silva_taxtab_spp.csv"), col.names=NA)

st12 <- readRDS("16S_seqtab.Rds")

# Construct Phylogenetic Tree
# Code sourced from: https://f1000research.com/articles/5-1492/v2

#Sort by higher-abundance taxa
st12 <- st12[,order(colSums(st12), decreasing=TRUE)]

library("DECIPHER")
# Multiple sequence alignment of sequence variants
seqs <- getSequences(st12)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
saveRDS(ma,"16S_alignment.Rds")
ma <- readRDS("16S_alignment.Rds")

# Build a maximum likelihood tree from a neighbor-joining tree
library("phangorn")
phang.align <- as.phyDat(as(alignment,"matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR,"16S_pml_tree.Rds")
saveRDS(fitGTR$tree,"16S_tree.Rds")

detach("package:phangorn", unload=TRUE)

