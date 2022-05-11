#This script is designed to get the ratios of repeats found over num proteins
#in the "query" (printed)
#Created By David E. Hufnagel on June 8 2012

import sys, os

rep = open(sys.argv[1]) #input 4col RM file
que = open(sys.argv[2]) #"query" fasta file





repSet = set()
for line in rep:
    lineLst = line.split("\t")
    repSet.add(lineLst[1])
    
repe = len(repSet)
tot = 0
for line in que:
    #print line
    if line.startswith(">"):
        #print line
        tot += 1
#print tot
#print repe
print float(repe) / float(tot)


rep.close()

