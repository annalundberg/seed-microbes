#! /usr/bin/env python3
# Filters OTU's from a mapping file, based on sequence headers in an otu fasta file

import sys

def main():
  map_file = open(sys.argv[1])
  fasta_file = open(sys.argv[2])
  filtered_map = open(sys.argv[3], 'w')

  select_otus = set()
  for line in fasta_file:
    if line[0] == '>':
      otu = line.rstrip()[1:]
      otu = otu.split(";size=")[0]
      select_otus.add(otu)
  
  fasta_file.close()
  
  for line in map_file:
     otu = line.split('\t')[0]
     if otu in select_otus:
       filtered_map.write(line)
  
  map_file.close()
  filtered_map.close()  

if __name__ == "__main__":
    main()   
