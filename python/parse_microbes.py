#! /usr/bin/env python3

### This script parses fastq files into 4 folders based on primer identity: 16S(799F-1193R), ITS1F/ITS2, shared ID's and rejects.####
### It reads the forward file only (usually better quality) and fills in the reverse file accordingly.
### I refined it using my data to maximize accepts but still distinguish between bacterial and fungal sequences.###
### This script also trims off variable length the ends of forward and reverse reads, so they start at the same universal sequence ###


# General usage: parse_microbes.py [forward.fastq] [OPTIONAL output diectory]
# With this script in the same folder as the fastqfiles, for example run 1, type: for f in *R1_run1.fastq; do python parse_microbes.py $f; done

import sys
import os


def main():
    fastqFname = str(sys.argv[1])
    fastqRname = fastqFname.replace('R1','R2') 
    writedir = "./"
    if len(sys.argv) > 2:
        writedir = str(sys.argv[2])
        if  writedir[-1] != '/':
            writedir=writedir+'/'
    if not os.path.exists(writedir):
        os.makedirs(writedir)
    F_16S = open(writedir+"16S_"+fastqFname, 'w')
    R_16S = open(writedir+"16S_"+fastqRname, 'w')
    F_ITS = open(writedir+"ITS_"+fastqFname, 'w')
    R_ITS = open(writedir+"ITS_"+fastqRname, 'w')
    if not os.path.exists(writedir+'assign_shared'):
        os.makedirs(writedir+'assign_shared')
    shared = open(writedir+"assign_shared/shared_"+fastqFname, 'w')
    if not os.path.exists(writedir+'assign_rejects/'):
        os.makedirs(writedir+'assign_rejects/')
    rejects = open(writedir+"assign_rejects/rejects.fastq, 'w')
    counter = 1
    bactF = 'AGATA'
    bactR = 'GTCA'
    fungF = 'GTAAA'
    fungR = 'GTTCTT'
    write16S = 0
    writeITS = 0
    reject = 0
    with open(fastqFname) as fastqF, open(fastqRname) as fastqR:
        for lineF,lineR in zip(fastqF, fastqR):
            if counter > 4:
                counter = 1
                reject = 0
            if lineF[0] == '@' and counter == 1:
                headF = lineF.rstrip()
                headR = lineR.rstrip()
                if headF[-2:-1] == '/': # remove bowtie2 header adjustment
                    headF =  headF[:-2]
                if headR[-2:-1] == '/': # remove bowtie2 header adjustment
                    headR =  headR[:-2] 
                if headF != headR:
                    sys.exit("F and R sequence headers do not match")
                head = headF+'\n'
            if counter == 2:
                seqF = lineF
                seqR = lineR
                if bactF in seqF[0:20] and bactR in seqR[0:16]:
                    write16S = 1
                    seqFstart = seqF.find(bactF, 0, 20)
                    seqRstart = seqR.find(bactR, 0, 16)
                    seqF = seqF[seqFstart:]
                    seqR = seqR[seqRstart:]
                elif fungF in seqF[0:20] and fungR in seqR[0:25]:
                    writeITS = 1
                    seqFstart = seqF.find(fungF, 0, 20)
                    seqRstart = seqR.find(fungR, 0, 25)
                    seqF = seqF[seqFstart:]
                    seqR = seqR[seqRstart:]
                else:
                    reject = 1
                if len(seqF) < 10 or len(seqR) < 10:
                    reject = 1
            if counter == 4:
                qualF = lineF
                qualR = lineR
                if reject == 0:
                    qualF = qualF[seqFstart:]
                    qualR = qualR[seqRstart:]
                    if write16S == 1 and writeITS ==1:
                        shared.write(head+seqF+'+\n'+qualF)
                        shared.write(head+seqR+'+\n'+qualR)
                        write16S = 0
                        writeITS = 0
                    elif write16S == 1:
                        F_16S.write(head+seqF+'+\n'+qualF)
                        R_16S.write(head+seqR+'+\n'+qualR)
                        write16S = 0
                    elif writeITS == 1:
                        F_ITS.write(head+seqF+'+\n'+qualF)
                        R_ITS.write(head+seqR+'+\n'+qualR)
                        writeITS = 0
                if reject == 1:
                    rejects.write(head+seqF+'+\n'+qualF)
                    rejects.write(head+seqR+'+\n'+qualR)
            counter = counter + 1 
    F_16S.close()
    R_16S.close()
    F_ITS.close()
    R_ITS.close()
    rejects.close()
    fastqF.close()
    fastqR.close()
    if 'shared' in locals() or 'shared' in globals():
        shared.close()
  

if __name__ == "__main__":
    main()        
