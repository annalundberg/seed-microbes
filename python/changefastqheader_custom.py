#! /usr/bin/env python3

import sys
import os

def main():
    '''Fxn obtains file using sys.argv index. Iterates through fastq file, each header is
    edited to fit custom format. Each fastq read with an edited header will be written into
    newfile in folder newheader.'''
    counter = 1
    with open(sys.argv[1], 'r') as fastqfile, open('newheader/'+fastqfile.name, 'w') as newfile:
        if not os.path.exists('newheader'):
            os.makedirs('newheader')
        for line in fastqfile:
            if line[0] == '@' and counter%4 == 1:
                line =  line[:-3]
                run = line.split("-")[1]
                line = line.split("-")[0]
                sample = line.split("_")[0]
                seq = line.split("_")[1]
                line = sample+'.'+run+'_'+seq+'\n'
            newfile.write(line)
            counter += 1
  
if __name__ == "__main__":
    main()        
