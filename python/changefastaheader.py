#! /usr/bin/env python3

import sys

def main():
    fastafile = open(sys.argv[1], 'r')
    newfile = list()
    for line in fastafile:
        if line[0] == '>':
            line =  line.rstrip()
            samplerun = line.split("_")[0]
            sample = samplerun.split(".")[0]
            run = samplerun.split(".")[1]
            read = line.split("_")[1][4:]
            line = sample+"_run"+run+".seqID="+read+'\n'
        newfile.append(line)    
    fastafile.close()           
    fastafile = open(sys.argv[1], 'w')
    for line in newfile:
        fastafile.write(line)    
       

if __name__ == "__main__":
    main()        
