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
#Run2
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
#This results in ~70% assembly or reads
system("pear -f bowtie2/microbes_R1.fastq -r bowtie2/microbes_R2.fastq -o pear_assembly/microbes -g 1 -j 8")

#Try m.a.p. minimum accepted probability for longer sequences
system("pear -f bowtie2/microbes_R1.fastq -r bowtie2/microbes_R2.fastq -o pear_assembly/microbes_map -g 2 -j 8 -v 30")



#################### Convert to Fasta and Qual format, remove bowtie header "/1' ########################

system('mothur "#fastq.info(fastq=pear_assembly/microbes.assembled.fastq);quit()"')
system("mkdir PDB_filtering")
system("mv pear_assembly/microbes.assembled.fasta PBD_filtering")
system("mv pear_assembly/microbes.assembled.qual PBD_filtering")


######################## DECIPHER de novo chimera detection
library(DECIPHER)


###################### Quality Filtering Step 1 ############################
#!!!THis script can use more automation, using system2() command and args=c() http://www.kdnuggets.com/2015/10/integrating-python-r-executing-part2.html
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/system2.html

setwd(paste(wd,'PBF_filtering',sep=''))
system("mkdir strict")
system("mkdir relaxed")

#The file is too large to filter on my computer - 15+ GB of RAM. I split the fasta and qual files into 4 equal parts based on line numbers. Still requires ~3.5 GB per run.
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
# Recombine parts and cleanup unnecessary output
system("rm microbes.assembled.part*")
setwd(paste(wd,'PBD_filtering/strict',sep=''))
system("cat *qc.good.fasta > microbes.hiquality.fasta")
system("cat *qc.good.qual > microbes.hiquality.qual")
system("cat *qc.bad.fasta > microbes.hqrejects.fasta")
system("rm microbes.assembled.part*")


setwd(paste(wd,'PBD_filtering/strict',sep=''))
system("cat *qc.good.fasta > microbes.hiquality.fasta")
system("cat *qc.bad.fasta > microbes.hqrejects.fasta")
system("rm microbes.assembled.part*")




### VSEARCH Pipeline
setwd(paste(wd,'vsearch',sep=''))

system("vsearch --threads 8 --fastq_eestats ../pear_assembly/microbes.assembled.fastq \ --output microbes.assembled.stats") 

# This is the current vsearch filter parameters, which is sent unto derep 
SWARM=$(which swarm)

#try with non-filtered seqs
system("vsearch --uchime_denovo microbes_assembled_dereplicated.fasta --threads 8 --sizeout --fasta_width 0 --nonchimeras microbes_nochimeras_denovo.fasta")

system("vsearch --uchime_ref microbes_nochimeras_denovo.fasta --db ../reference_database/uchime_16S_ITS.fasta --sizeout --fasta_width 0 --nonchimeras microbes_nochimeras_denovo_ref.fasta --threads 8")



# Identify Chimeras
system("identify_chimeric_seqs.py -m usearch61 -r ../reference_database/16S_ITS_97_otus.fasta -i strict/microbes.hiquality.fasta -o ../usearch61_chimera_check/")
system("identify_chimeric_seqs.py -m usearch61 -r ../reference_database/16S_ITS_97_otus.fasta -i relaxed/microbes.quality.fasta -o ../usearch_otus/usearch61_chimera_check_relaxed/")

setwd(wd)
system("mkdir usearch_otus")
setwd(paste(wd,'usearch_otus',sep=''))

system("filter_fasta.py -f ../PBD_filtering/strict/microbes.hiquality.fasta -o microbes_strict_nochime.fasta -s usearch61_chimera_check/non_chimeras.txt")

#Start with high-quality reads to establish OTUS using de novo approach -k suppresses de novo chimera detection, -x suppresses reference chimera detection
system('pick_otus.py -i microbes_strict_nochime.fasta -m usearch61 -o usearch_denovo_otus/')
system('pick_otus.py -i ../PBD_filtering/strict/microbes.hiquality.fasta -m usearch -o usearch_denovo_otus/ --db_filepath ../reference_database/16S_ITS_97_otus.fasta')

system('pick_otus.py -i ../PBD_filtering/relaxed/microbes.quality.fasta -m usearch_ref -r usearch_denovo_otus/enumerated_otus.fasta -o usearch_open_ref_otus/ --db_filepath ../reference_database/16S_ITS_97_otus.fasta')

system("pick_open_reference_otus.py -i ../PBF_filtering/relaxed/microbes.quality.fasta -m usearch61 -r usearch_denovo_otus/enumerated_otus.fasta -o usearch_open_ref_otus/ --parallel --suppress_taxonomy_assignment --suppress_align_and_tree")
system("identify_chimeric_seqs.py -m usearch61 -r ../reference_database/16S_ITS_97_otus.fasta -i usearch_open_ref_otus/rep_set.fna -o usearch_open_ref_otus/")



