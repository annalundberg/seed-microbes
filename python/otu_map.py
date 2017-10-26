#! /usr/bin/env python3
# Takes a .uc file from clustered sequences, and a .uc file from dereplicated sequences
# Outputs a map file, with all sequences listed for each cluster

import sys

def main():
  uc_derep = open(sys.argv[1])
  uc_clust = open(sys.argv[2])
  map_file = open(sys.argv[3], 'w')

  # Create a dictionary of the clusters, with representative seeds of each cluster as each entry reference
  # add all cluster members as values of the reference
  
  cdict = {}
  for line in uc_clust:
    cluster = int(line.split("\t")[1])
    query = line.split("\t")[8]
    query = query.split(";size=")[0]
    target = target = line.split("\t")[9]
    target = target.split(";size=")[0]
    if line.split("\t")[0] == 'S':
      cdict[query] = [cluster,query]
    if line.split("\t")[0] == 'H':
      cdict[target].append(query)
  
  uc_clust.close()
  
  # Create dictionary of derep, where representative sequence is the key, and mathcing sequences are the values
  
  ddict = {}
  for line in uc_derep:
    cluster = int(line.split("\t")[1])
    query = line.split("\t")[8]
    target = (line.split("\t")[9]).rstrip()
    if line.split("\t")[0] == 'S':
      ddict[query] = [query]
    if line.split("\t")[0] == 'H':
      ddict[target].append(query) 

  uc_derep.close()
  
  # Parse through the clustering dictionary and add all dereplicated sequences to each cluster
  
  # for otu in cdict:
  #   for derep in cdict[otu]:
  #     if isinstance(derep,int):
  #       otuID = derep + 1
  #       otuID = 'OTU_'+str(otuID)
  #       line = otuID+'\t'
  #     else:
  #       reps = ddict[derep]
  #       for rep in reps:
  #         line = line+rep+'\t'
  #   map_file.write(line+'\n')
    
  otu_list = []
  for num in range(len(cdict)):
    otu_list.append([num])
    
  for ref in cdict:
    otuID = cdict[ref][0]
    for derep in cdict[ref][1:len(cdict[ref])]:
      reps = ddict[derep]
      for rep in reps:
        otu_list[otuID].append(rep)

  for otu in otu_list:
    otuID = 'OTU_'+str(otu[0]+1)
    line = otuID+'\t'
    for seq in otu[1:len(otu)]:
      line = line + seq + '\t'
    map_file.write(line+'\n')  
      
  map_file.close()
    

if __name__ == "__main__":
    main()   
