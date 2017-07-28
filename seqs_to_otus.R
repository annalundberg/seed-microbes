library(rPython)

#Set working directory
wd <- "/home/lucas/Main/Illumina/pipeline/"
#Set wd to location of this script
#wd <- dirname(parent.frame(2)$ofile)
#setwd(wd)

setwd(wd)
#Add the python directory to the system PATH variable so it can run custom python scripts
ppath <- paste(wd,'python',sep="")
syspath <- paste(Sys.getenv("PATH"), ppath, sep=":")
Sys.setenv(PATH=syspath)
system("echo 'New folder added to PATH variable:' $PATH")

#Remultiplex the sequences with custom Python script: remultiplex.py <input files> <output file> <run#>
# The script relabels the header of each sequence: <Sample#>_<Sequence#>/<Run#>

setwd(paste(wd,"run1/run1_raw_fastq", sep=""))
system('for f in *R1_001.fastq; do remultiplexfastq.py $f R1_run1.fastq 1; done')
system('for f in *R2_001.fastq; do remultiplexfastq.py $f R2_run1.fastq 1; done')
setwd(paste(wd,"run2/run2_raw_fastq", sep=""))
system('for f in *R1_001.fastq; do remultiplexfastq.py $f R1_run2.fastq 2; done')
system('for f in *R2_001.fastq; do remultiplexfastq.py $f R2_run2.fastq 2; done')

#Then concatenate the combined runs into a single data file
setwd(wd)
system('cat run1/run1_raw_fastq/R1_run1.fastq run2/run2_raw_fastq/R1_run2.fastq > complete_R1.fastq')
system('cat run1/run1_raw_fastq/R2_run1.fastq run2/run2_raw_fastq/R2_run2.fastq > complete_R2.fastq')

############### Use Bowtie2 to remove host Maize DNA ###########################
setwd(paste(wd,'bowtie',sep=''))

# 1) download Zea mays genome and create a new index called "maize"
system("bowtie2-build Zea_mays.AGPv4.dna.chromosome.1.fa maize")
# 2) bowtie2 mapping against host hg19, keep both mapped and unmapped reads (paired-end reads)
system("bowtie2 -x maize -1 ../complete_R1.fastq -2 ../complete_R2.fastq -S complete_mapped_and_unmapped.sam")
# 3) convert file .sam to .bam
system("samtools view -bS complete_mapped_and_unmapped.sam > complete_mapped_and_unmapped.bam")
# 4) filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
system("samtools view -b -f 12 -F 256 complete_mapped_and_unmapped.bam > complete_bothEndsUnmapped.bam")
# -f 12     Extract only (-f) alignments with both reads unmapped: <read unmapped><mate unmapped>
# -F 256   Do not(-F) extract alignments which are: <not primary alignment>
# 5)  split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
system("samtools sort -n complete_bothEndsUnmapped.bam complete_bothEndsUnmapped_sorted.bam")
system("bedtools bamtofastq -i SAMPLE_bothEndsUnmapped_sorted.bam -fq complete_nonhost_r1.fastq -fq2 complete_nonhost_r2.fastq")

##############
setwd(wd)
system("mkdir pear")
system("pear -f SAMPLE_r1.fastq -r SAMPLE_r2.fastq -o pear/ -g 1 -j 8")
system("pear -f SAMPLE_r1.fastq -r SAMPLE_r2.fastq -o pear2 -g 2 -j 8 -p 1.0 -b 24")

