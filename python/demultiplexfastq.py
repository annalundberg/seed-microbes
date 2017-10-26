#! /usr/bin/env python3
# Takes a .fastq file and reads headers for sample names, and generates separate .fastq files for each sample
# Assumes that file names have the following format:[AMPLICON]_[Sample]_[run#].fastq
# Required input fastq header format: Sample[SampleName]."[Run#]_Sequence[Sequence#]

import sys
import os

def main():
    counter = 1
    prevsample = ''
    fastqfile = open(sys.argv[1], 'r')
    outdir = "./"
    if len(sys.argv) > 2:
      outdir = sys.argv[2]
      if outdir[-1] != '/':
        outdir = outdir + '/'
    amplicon = fastqfile.name[0:3]
    if not os.path.exists(outdir+'combined/'+amplicon):
        os.makedirs(outdir+'combined/'+amplicon)
    if not os.path.exists(outdir+'run1/'+amplicon):
        os.makedirs(outdir+'run1/'+amplicon)
    if not os.path.exists(outdir+'run2/'+amplicon):
        os.makedirs(outdir+'run2/'+amplicon)
    for line in fastqfile:
        if counter > 4:
            counter = 1
        if line[0] == '@' and counter == 1:
            samplerun = line.split("_")[0][1:]
            sample = samplerun.split(".")[0]
            run = samplerun.split(".")[1]
        if samplerun != prevsample:
            prevsample = samplerun
            fastqsample = open(outdir+'combined/'+amplicon+"/"+amplicon+"_"+str(sample)[6:]+'_combined.fastq', 'a')
            if run=="1":
                fastqsample1 = open(outdir+'run1/'+amplicon+"/"+amplicon+"_"+str(sample)[6:]+'_run1.fastq', 'a')
            if run=="2":
                fastqsample2 = open(outdir+'run2/'+amplicon+"/"+amplicon+"_"+str(sample)[6:]+'_run2.fastq', 'a')
        fastqsample.write(line)
        if run=="1":
            fastqsample1.write(line)
        if run=="2":
            fastqsample2.write(line)
        counter = counter + 1  
    fastqsample.close()
    fastqsample1.close()
    fastqsample2.close()
    fastqfile.close()       
        

if __name__ == "__main__":
    main()        
