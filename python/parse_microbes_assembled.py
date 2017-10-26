#! /usr/bin/env python3

### This script parses fastq files into 4 folders based on primer identity: 16S(799F-1193R), ITS1, shared ID's and rejects.####
# It reads the forward file only (usually better quality) and fills in the reverse file accordingly.
### I refined it using my data to maximize accepts but still distinguish between bacterial and fungal sequences.###

# With this script in the same folder as the fastqfiles, for example run 1, type: for f in *R1_run1.fastq; do python parse_microbes.py $f; done

import sys
import os
import linecache

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def main():
    fastqF = open(sys.argv[1], 'r')
    if not os.path.exists('16S'):
        os.makedirs('16S')
    if not os.path.exists('ITS'):
        os.makedirs('ITS')
    R1_16S = open("16S/16S_"+fastqF.name, 'w')
    R1_ITS = open("ITS/ITS_"+fastqF.name, 'w')
    counter = 1
    linecount=1
    bact = 'AGATAC'
    fung = 'GTAAA'
    write16S = 0
    writeITS = 0
    seq=""
    qual=""
    for line in fastqF:
        if counter > 4:
            counter = 1
        if line[0] == '@' and counter == 1:
            head = line
        if counter == 2:
            seqF = line
            if bact in seqF[2:20]:
                write16S = 1
            if fung in seqF[2:20]:
                writeITS = 1
        if counter == 4:
            qualF = line
            if write16S == 1 and writeITS ==1:
                if not os.path.exists('assign_shared'):
                    os.makedirs('assign_shared')
                sharedF = open("assign_shared/shared_"+fastqF.name, 'a')
                sharedF.write(head+seqF+'+\n'+qualF)
                write16S = 0
                writeITS = 0
            if write16S == 1:
                R1_16S.write(head+seqF+'+\n'+qualF)
                write16S = 0
            elif writeITS == 1:
                R1_ITS.write(head+seqF+'+\n'+qualF)
                writeITS = 0
            else:
                if not os.path.exists('assign_rejects'):
                    os.makedirs('assign_rejects')
                rejectsF = open("assign_rejects/rejects_"+fastqF.name, 'a')
                rejectsF.write(head+seqF+'+\n'+qualF)
        counter = counter + 1 
        linecount = linecount + 1 
    fastqF.close() 
    R1_16S.close()
    R1_ITS.close()
    if 'rejectsF' in locals() or 'rejectsF' in globals():
        rejectsF.close()
    if 'sharedF' in locals() or 'sharedF' in globals():
        sharedF.close()

if __name__ == "__main__":
    main()        
