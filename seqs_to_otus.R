library(rPython)

#Set working directory

wd <- "/home/lucas/Main/Illumina/seed-microbes/"
setwd(wd)

#system("mkdir python")
#system("for f in *.py; do mv $f python; done")
#system("mkdir R")


#Add the python directory to the system PATH variable so it can run custom python scripts
ppath <- paste(wd,'python',sep="")
syspath <- paste(Sys.getenv("PATH"), ppath, sep=":")
Sys.setenv(PATH=syspath)
system("echo 'New folder added to PATH variable:' $PATH")

#Remultiplex the sequences with custom Python script: remultiplex.py <input files> <output file> <run#>
# The script relabels the header of each sequence: <Sample#>_<Sequence#>/<Run#>
#Run1
system("mkdir run1")
system("mv seed-microbes_run1_raw_fastq.tar.gz run1")
setwd(paste(wd,'run1',sep=""))
system("tar -xzf seed-microbes_run1_raw_fastq.tar.gz")
setwd(paste(wd,"run1", sep=""))
system('for f in *R1_001.fastq; do remultiplexfastq.py $f run1_R1.fastq 1; rm $f; done')
system('for f in *R2_001.fastq; do remultiplexfastq.py $f run1_R2.fastq 1; rm $f; done')
#Run1
system("mkdir run2")
system("mv seed-microbes_run2_raw_fastq.tar.gz run2")
setwd(paste(wd,'run2',sep=""))
system("tar -xzf seed-microbes_run2_raw_fastq.tar.gz")system('for f in *R2_001.fastq; do remultiplexfastq.py $f run1_R2.fastq 1; done')
setwd(paste(wd,"run2", sep=""))
system('for f in *R1_001.fastq; do remultiplexfastq.py $f run2_R1.fastq 2; rm $f; done')
system('for f in *R2_001.fastq; do remultiplexfastq.py $f run2_R2.fastq 2; rm $f; done')

#Then concatenate the combined runs into a single data file
setwd(wd)
system('cat run1/run1_R1.fastq run2/run2_R1.fastq > complete_R1.fastq')
system('cat run1/run1_R2.fastq run2/run2_R2.fastq > complete_R2.fastq')

############### Use Bowtie2 to remove host Maize DNA ###########################
setwd(paste(wd,'bowtie2',sep=''))

# 1) download Zea mays genome and create a new index called "maize"
system("bowtie2-build Zea_mays.AGPv4.dna.chromosome.1.fa maize")
# 2) bowtie2 mapping against host hg19, keep both mapped and unmapped reads (paired-end reads)
#This step took 7.5 hours, but without the multithread option
system("bowtie2 -x maize -1 ../complete_R1.fastq -2 ../complete_R2.fastq -S complete_mapped_and_unmapped.sam -p 8")
# 3) convert file .sam to .bam
#-@8 signifies 8 threads
system("samtools view -bS -@ 8 complete_mapped_and_unmapped.sam > complete_mapped_and_unmapped.bam")
# 4) filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
system("samtools view -b -f 12 -F 256 -@ 8 complete_mapped_and_unmapped.bam > complete_bothEndsUnmapped.bam")
# -f 12     Extract only (-f) alignments with both reads unmapped: <read unmapped><mate unmapped>
# -F 256   Do not(-F) extract alignments which are: <not primary alignment>
# 5)  split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
system("samtools sort -n complete_bothEndsUnmapped.bam complete_bothEndsUnmapped_sorted")
system("bedtools bamtofastq -i complete_bothEndsUnmapped_sorted.bam -fq microbes_R1.fastq -fq2 microbes_R2.fastq")

system("changefastqheader.py microbes_R1.fastq microbes_R1_newheader.fastq")
system("mv microbes_R1_newheader.fastq microbes_R1.fastq")
system("changefastqheader.py microbes_R2.fastq microbes_R2_newheader.fastq")
system("mv microbes_R2_newheader.fastq microbes_R2.fastq")

############## Assemble paired reads using PEAR ############################
setwd(wd)
system("mkdir pear_assembly")
system("pear -f bowtie2/microbes_R1.fastq -r bowtie2/microbes_R2.fastq -o pear_assembly/microbes -g 1 -j 8")


#################### Convert to Fasta and Qual format, remove bowtie header "/1' ########################

system('mothur "#fastq.info(fastq=pear_assembly/microbes.assembled.fastq);quit()"')
system("mkdir PDB_filtering")
system("mv pear_assembly/microbes.assembled.fasta PBD_filtering")
system("mv pear_assembly/microbes.assembled.qual PBD_filtering")



###################### Quality Filtering Step 1 ############################
#!!!THis script can use more automation, using system2() command and args=c() http://www.kdnuggets.com/2015/10/integrating-python-r-executing-part2.html
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/system2.html

setwd(paste(wd,'PBD_filtering',sep=''))
system("mkdir strict")
system("mkdir relaxed")

#The file is too large to filter - 15+ GB of RAM. I split the fasta and qual files into 4 equal parts based on line numbers. Still requires ~3.5 GB per run.
system("split -l$((`wc -l < microbes.assembled.fasta`/4)) microbes.assembled.fasta microbes.assembled.part --additional-suffix='.fasta' -da 1")
system("split -l$((`wc -l < microbes.assembled.qual`/4)) microbes.assembled.qual microbes.assembled.part --additional-suffix='.qual' -da 1")

# Use default setings for Poisson Binomial Distribution filtering to acquire highest quality reads for OTU picking
system('moira.py --forward_fasta=microbes.assembled.part0.fasta --forward_qual=microbes.assembled.part0.qual --processors 8 --output_prefix=strict/microbes.assembled.part0')
system('moira.py --forward_fasta=microbes.assembled.part1.fasta --forward_qual=microbes.assembled.part1.qual --processors 8 --output_prefix=strict/microbes.assembled.part1')
system('moira.py --forward_fasta=microbes.assembled.part2.fasta --forward_qual=microbes.assembled.part2.qual --processors 8 --output_prefix=strict/microbes.assembled.part2')
system('moira.py --forward_fasta=microbes.assembled.part3.fasta --forward_qual=microbes.assembled.part3.qual --processors 8 --output_prefix=strict/microbes.assembled.part3')
# Relax the quality filtering settings to keep more sequences
system('moira.py --forward_fasta=microbes.assembled.part0.fasta --forward_qual=microbes.assembled.part0.qual --processors 8 --output_prefix=relaxed/microbes.assembled.part0 --uncert=0.02')
system('moira.py --forward_fasta=microbes.assembled.part1.fasta --forward_qual=microbes.assembled.part1.qual --processors 8 --output_prefix=relaxed/microbes.assembled.part1 --uncert=0.02')
system('moira.py --forward_fasta=microbes.assembled.part2.fasta --forward_qual=microbes.assembled.part2.qual --processors 8 --output_prefix=relaxed/microbes.assembled.part2 --uncert=0.02')
system('moira.py --forward_fasta=microbes.assembled.part3.fasta --forward_qual=microbes.assembled.part3.qual --processors 8 --output_prefix=relaxed/microbes.assembled.part3 --uncert=0.02')


setwd(paste(wd,'PBD_filtering/relaxed',sep=''))
system("cat *qc.good.fasta > microbes.quality.fasta")
system("cat *qc.good.qual > microbes.quality.qual")
system("cat *qc.bad.fasta > microbes.rejects.fasta")
system("rm microbes.assembled.part*")

# Identify Chimeras
system("identify_chimeric_seqs.py -m usearch61 -r ../../reference_database/16S_ITS_97_otus.fasta -i /strict/microbes.quality.fasta -o ../usearch61_chimera_check/")



system('moira.py --forward_fasta=relaxed/microbes.assembled.part0.fasta --forward_qual=microbes.assembled.part0.qual --processors 8 --output_prefix=strict2/microbes.assembled.part0')
system('moira.py --forward_fasta=relaxed/microbes.assembled.part1.fasta --forward_qual=microbes.assembled.part1.qual --processors 8 --output_prefix=strict2/microbes.assembled.part1')
system('moira.py --forward_fasta=relaxed/microbes.assembled.part2.fasta --forward_qual=microbes.assembled.part2.qual --processors 8 --output_prefix=strict2/microbes.assembled.part2')
system('moira.py --forward_fasta=relaxed/microbes.assembled.part3.fasta --forward_qual=microbes.assembled.part3.qual --processors 8 --output_prefix=strict2/microbes.assembled.part3')



# Recombine parts and cleanup unnecessary output
system("rm microbes.assembled.part*")
setwd(paste(wd,'PBD_filtering/strict',sep=''))
system("cat *qc.good.fasta > microbes.hiquality.fasta")
system("cat *qc.bad.fasta > microbes.hqrejects.fasta")
system("rm microbes.assembled.part*")



setwd(wd)
system("mkdir' otus")
setwd(paste(wd,'otus',sep=''))


#Start with high-quality reads to establish OTUS using de novo approach -k suppresses de novo chimera detection, -x suppresses reference chimera detection
system('pick_otus.py -i ../PBD_filtering/quality/microbes.hiquality.fasta -m usearch -o usearch_denovo_otus/ -k -x')



system('pick_otus.py -i ../PBD_filtering/quality/microbes.hiquality.fasta -m usearch -o usearch_denovo_otus/ -k -x')


