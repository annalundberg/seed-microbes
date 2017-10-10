#! /usr/bin/env python3


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
    filtF = open(writedir+"eq_"+fastqFname, 'w')
    filtR = open(writedir+"eq_"+fastqRname, 'w')
    
    search = 0

    linecount = 1
    counter = 1
    with open(fastqFname) as fastqF, open(fastqRname) as fastqR:
        for lineF,lineR in zip(fastqF, fastqR):
            if counter > 4:
              counter = 1
            if counter == 1:
              if lineF != line R:
                search = 1
                qF = lineF
                qR = lineR
                saveF <- lineF
                saveR <- lineR
                continue
            if search == 0:
              filtF.write(line)
              filtR.write(line)
            if search == 1:
              if counter != 1
                saveF <- saveF + lineF
                saveR <- saveR + lineR
              if counter == 1
                if qF == :
                if lineR in saveR:
              
              
            counter = counter + 1 
            linecount = linecount + 1     
    R1_16S.close()
    R2_16S.close()
    R1_ITS.close()
    R2_ITS.close()
    rejectsF.close()
    rejectsR.close()        
    fastqF.close()   
    if 'sharedF' in locals() or 'sharedF' in globals():
        sharedF.close()
        sharedR.close()    

if __name__ == "__main__":
    main()        
