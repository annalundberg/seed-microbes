#!/bin/sh


THREADS=8
REF=../reference_database/16S_ITS_97_otus.fasta
PERL=$(which perl)
VSEARCH=$(which vsearch)
WD=/home/lucas/Main/Illumina/seed-microbes
RUNS=2
PATH=$PATH:$WD/python

#Remultiplex the sequences with custom Python script: remultiplex.py <input files> <output file> <run#>
# The script relabels the header of each sequence: "Sample"[Sample#].[Run#]_"Sequence["Sequence#]


for i in $(seq 1 $RUNS)
do     
    cd $WD/raw_fastq/run$i
    for f in *R1_001.fastq.gz    
    do 
        gzip -d $f
        remultiplexfastq.py ${f%.*} ../prefilter/run"$i"_R1.fastq $i
    done 
    for f in *R2_001.fastq.gz
    do 
        gzip -d $f
        remultiplexfastq.py ${f%.*} ../prefilter/run"$i"_R2.fastq $i
    done
done 

#Then concatenate the combined runs into a single data file
cd $WD/prefilter/

cat run1_R1.fastq run2_R1.fastq > complete_R1.fastq
cat run1_R2.fastq run2_R2.fastq > complete_R2.fastq

##### Use Bowtie2 to remove host Maize DNA  ####
cd $WD/prefilter/bowtie2

## 1) download Zea mays genome and create a new index called "maize"
#bowtie2-build Zea_mays.AGPv4.dna.chromosome.1.fa maize
# 2) bowtie2 mapping against host hg19, keep both mapped and unmapped reads (paired-end reads)
# This step took 7.5 hours, but without the multithread option
bowtie2 -x maize -1 ../complete_R1.fastq -2 ../complete_R2.fastq -S complete_mapped_and_unmapped.sam -p $THREADS
# 3) convert file .sam to .bam
#-@8 signifies 8 threads
samtools view -bS -@ 8 complete_mapped_and_unmapped.sam > microbes_mapped_and_unmapped.bam
# 4) filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -b -f 12 -F 256 -@ 8 microbes_mapped_and_unmapped.bam > microbes_bothEndsUnmapped.bam
# -f 12     Extract only (-f) alignments with both reads unmapped: <read unmapped><mate unmapped>
# -F 256   Do not(-F) extract alignments which are: <not primary alignment>
# 5)  split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -n microbes_bothEndsUnmapped.bam microbes_bothEndsUnmapped_sorted
bedtools bamtofastq -i microbes_bothEndsUnmapped_sorted.bam -fq ../microbes_R1.fastq -fq2 ../microbes_R2.fastq

cd $WD/prefilter/

# Seperate 16S from ITS sequences using forward read .fastq files
# General usage: parse_microbes.py [forward.fastq] [OPTIONAL output diectory]
parse_microbes.py microbes_R1.fastq

# # Trim off sequences where quality drops off, before PEAR assembly
# **** I found that this does not improve PEAR alignment in my sequences. Will just quality-filter post-alignment ****

# $VSEARCH --fastq_filter 16S_microbes_R1.fastq \
#     --threads $THREADS \
#     --fastq_stripleft 10 \
#     --fastq_minlen 100 \
#     --fastq_truncqual 2 \
#     --fastqout 16S_microbes_R1_filtered.fastq \
#     --fasta_width 0 \
#     --fastaout_discarded 16S_microbes_R1_discarded.fasta
#     
#     
# $VSEARCH --fastq_filter 16S_microbes_R2.fastq \
#     --threads $THREADS \
#     --fastq_stripleft 10 \
#     --fastq_minlen 100 \
#     --fastq_truncqual 2 \
#     --fastqout 16S_microbes_R2_filtered.fastq \
#     --fasta_width 0 \
#     --fastaout_discarded 16S_microbes_R2_discarded.fasta
# 
# # replace empty lines in .fasta files with "N"
# sed -i -e 's/^$/N/' 16S_microbes_R1_discarded.fasta
# sed -i -e 's/^$/N/' 16S_microbes_R2_discarded.fasta
# 
# filter_fasta.py -f 16S_microbes_R1_filtered.fastq \
#   -o 16S_microbes_R1_filtered2.fastq \
#   -a 16S_microbes_R2_discarded.fasta -n
# filter_fasta.py -f 16S_microbes_R2_filtered.fastq \
#   -o 16S_microbes_R2_filtered2.fastq \
#   -a 16S_microbes_R1_discarded.fasta -n


############## Assemble paired reads using PEAR ############################ 
cd $WD/pear_assembly

pear -f ../prefilter/16S_microbes_R1.fastq -r ../prefilter/16S_microbes_R2.fastq -o ./16S -g 1 -j 8 -n 380 -m 460
# Assembled reads ...................: 2,716,456 / 4,113,314 (66.041%)
# Discarded reads ...................: 99 / 4,113,314 (0.002%)
# Not assembled reads ...............: 1,396,759 / 4,113,314 (33.957%)

# With stricter assembly size and trimming the spacers at the end ofeach sequence
# Assembled reads ...................: 1,790,811 / 3,288,309 (54.460%)
# Discarded reads ...................: 0 / 3,288,309 (0.000%)
# Not assembled reads ...............: 1,497,498 / 3,288,309 (45.540%)

pear -f ../prefilter/ITS_microbes_R1.fastq -r ../prefilter/ITS_microbes_R2.fastq -o ./ITS -g 1 -j 8 -n 50 -m 600
# Assembled reads ...................: 678,230 / 708,041 (95.790%)
# Discarded reads ...................: 34 / 708,041 (0.005%)
# Not assembled reads ...............: 29,777 / 708,041 (4.206%)

# Assembled reads ...................: 666,831 / 696,127 (95.792%)
# Discarded reads ...................: 0 / 696,127 (0.000%)
# Not assembled reads ...............: 29,296 / 696,127 (4.208%)

### Here, we diverge between dada2 and other clustering methods


# Demultiplex for dada2
ls
demultiplexfastq.py 16S.assembled.fastq ../dada2
demultiplexfastq.py ITS.assembled.fastq ../dada2


# Run DADA2 R script to
# 1. Filter assembled sequences with EE (expected error <= 2), OUTPUT amplicon_samplefastq.gz
# 2. Estimate error rates, OUTPUT: err.Rds
# 3. Learn sequence variants 


######  VSEARCH clustering ############
THREADS=8
REF=../reference_database/UCHIME/16S_UCHIME_gold.fa
VSEARCH=$(which vsearch)
WD=/home/lucas/Main/Illumina/seed-microbes
RUNS=2
PATH=$PATH:$WD/python
cd ~/Main/Illumina/seed-microbes/vsearch

### Filter and Dereplicate

$VSEARCH --fastq_filter ../pear_assembly/16S.assembled.fastq \
    --threads $THREADS \
    --fastq_minlen 380 \
    --fastq_maxlen 460 \
    --fastq_maxee 2 \
    --fastq_maxns 0 \
    --fastaout 16S.assembled.filtered.fasta \
    --fasta_width 0 \
    --fastaout_discarded 16S.assembled.discarded.fasta

echo
echo Dereplicate across samples

$VSEARCH --threads $THREADS \
    --derep_fulllength 16S.assembled.filtered.fasta \
    --sizeout \
    --fasta_width 0 \
    --uc 16S.derep.uc \
    --output 16S.derep.fasta


####### VSEARCH Cluster into OTUs before chimera detection ##########


$VSEARCH --threads $THREADS \
    --cluster_size 16S.derep.fasta \
    --id 0.97 \
    --minuniquesize 2 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc 16S.all.clustered.uc \
    --relabel OTU_ \
    --centroids 16S.all.otus.fasta
#    --otutabout all.otutab.txt


# Create OTU mapping file using custom script using uclust files generated from clustering and dereplication
otu_map.py 16S.derep.uc 16S.all.clustered.uc 16S.all.otu.map.txt

## Refine OTUs

# remove singleton otus
$VSEARCH --sortbysize 16S.all.otus.fasta --output 16S.otus.fasta --minsize 2

# Remove chimeric OTUs

$VSEARCH --threads $THREADS \
    --uchime_denovo 16S.otus.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --chimeras 16S.otus.denovo.chimeras.fasta \
    --nonchimeras 16S.otus.denovo.nc.fasta

echo OTUs after denovo chimera detection: $(grep -c "^>" 16S.otus.denovo.nc.fasta)
# OTUs after denovo chimera detection: 15587

$VSEARCH --threads $THREADS \
    --uchime_ref 16S.otus.denovo.nc.fasta \
    --db $REF \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --chimeras 16S.otus.denovo.ref.chimeras.fasta \
    --nonchimeras 16S.otus.denovo.ref.nc.fasta

# Found 787 (5.0%) chimeras, 14739 (94.6%) non-chimeras,
# and 61 (0.4%) borderline sequences in 15587 unique sequences.
# Taking abundance information into account, this corresponds to
# 29253 (2.6%) chimeras, 1014111 (88.6%) non-chimeras,
# and 101146 (8.8%) borderline sequences in 1144510 total sequences.

echo OTUs after reference-based chimera detection: $(grep -c "^>" 16S.otus.denovo.ref.nc.fasta)
# OTUs after reference-based chimera detection: 14739

# 
filter_otumap.py 16S.all.otu.map.txt 16S.otus.denovo.ref.nc.fasta 16S.final.otu.map.txt


# Create OTU table
make_otu_table.py -i 16S.final.otu.map.txt -o 16S.otu.table.biom


###### VSEARCH chimera detection before OTU generation ###########
cd ~/Main/Illumina/seed-microbes/vsearch/prefilter_chimeras


echo
echo Precluster at 98% before chimera detection

$VSEARCH --threads $THREADS \
    --cluster_size ../16S.derep.fasta \
    --id 0.98 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc 16S.preclustered.uc \
    --centroids 16S.preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>" all.preclustered.fasta)

echo
echo De novo chimera detection

$VSEARCH --threads $THREADS \
    --uchime_denovo 16S.preclustered.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras 16S.preclustered.denovo.nonchimeras.fasta

# Found 126963 (33.2%) chimeras, 255069 (66.6%) non-chimeras,
# and 761 (0.2%) borderline sequences in 382793 unique sequences.
# Taking abundance information into account, this corresponds to
# 209561 (14.7%) chimeras, 1213244 (85.2%) non-chimeras,
# and 1888 (0.1%) borderline sequences in 1424693 total sequences

echo Unique sequences after denovo chimera detection: $(grep -c "^>" 16S.preclustered.denovo.nonchimeras.fasta)
# Unique sequences after denovo chimera detection: 255069

$VSEARCH --threads $THREADS \
    --uchime_ref 16S.preclustered.denovo.nonchimeras.fasta \
    --db $REF \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --chimeras 16S.preclustered.ref.chimeras.fasta \
    --nonchimeras 16S.preclustered.ref.nonchimeras.fasta

# Found 11811 (4.6%) chimeras, 242489 (95.1%) non-chimeras,
# and 769 (0.3%) borderline sequences in 255069 unique sequences.
# Taking abundance information into account, this corresponds to
# 33791 (2.8%) chimeras, 1091549 (90.0%) non-chimeras,
# and 87904 (7.2%) borderline sequences in 1213244 total sequences.

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" 16S.preclustered.ref.nonchimeras.fasta)
# Unique sequences after reference-based chimera detection: 242489

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

uncluster.py ../16S.derep.fasta 16S.preclustered.uc 16S.preclustered.ref.nonchimeras.fasta 16S.derep.nonchimeras.fasta

echo
echo Cluster at 97% and relabel with OTU_n, generate OTU table

$VSEARCH --threads $THREADS \
    --cluster_size 16S.derep.nonchimeras.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel OTU_ \
    --centroids 16S.otus.fasta 
#    --otutabout all.otutab.txt
    
# This script only uses abundance 



echo
echo Number of OTUs: $(grep -c "^>" 16S.otus.fasta)
# Number of OTUs: 98869

# In 918207 sequences

