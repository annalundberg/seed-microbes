#! /usr/bin/env python3
# Uses the following arguements:
# 1. original fasta file before clustering
# 2. .uc file from clustering
# 3. fasta file after clustering

import sys

def main():
  orig_fasta = open(sys.argv[1])
  uclust = open(sys.argv[2])
  select_fasta = open(sys.argv[3])
  select_headers = set()
  out_fasta = open(sys.argv[4], 'w')

  for line in select_fasta:
    if line[0] == '>':
      head = line.rstrip()[1:]
      head = head.split(";size=")[0]
      select_headers.add(head)
  
  select_fasta.close()
  
  for line in uclust:
    if line.split("\t")[0] == 'H':
      seed = line.split("\t")[9]
      seed = seed.split(";size=")[0]
      if seed in select_headers:
        newseq = line.split("\t")[8]
        newseq = newseq.split(";size=")[0]
        select_headers.add(newseq)
  
  uclust.close()
  
  write = 0
  query = ''
  for line in orig_fasta:
    if line[0] == '>':
      head = line.rstrip()[1:]
      query = head.split(";size=")[0]
      if query in select_headers:
        write = 1  
        continue
    if write == 1:
      out_fasta.write(">"+head+'\n'+line)
      write = 0
  
  orig_fasta.close()
  out_fasta.close()

if __name__ == "__main__":
    main()   
  
