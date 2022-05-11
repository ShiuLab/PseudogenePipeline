#This script is intended to take an ortholog file with just two gene names
#and remove lines from a blast to create an m8 blast results file with only
#the orthologous gene pairs specified by the ortholog file
#Created by David E. Hufnagel on 3-8-2012

"""
Algorithm:
    1) iterate through ortholog file
        a) iterate through the blast file to search for the same pair
        b) when found output the blast line
        
"""

import sys

blast = open(sys.argv[1])     # the full m8 blast output file
ortho = open(sys.argv[2])     # the orthologs file
out = open(sys.argv[3], "w")  # the m8 blast output file with only ortologs

#put the blast file into a tuple of tuples
blastLst = []
for l in blast:
    blastLst.append(l)
    
#get and write the orthologous blast lines
oCnt = 0
for line in ortho:
    lineLst = line.split("\t")
    oName1 = lineLst[0] #gene names 1 and 2 from ortholog file
    oName2 = lineLst[1][:-1]
    for ln in blastLst:
        lnLst = ln.split("\t")
        bName1 = lnLst[0] #gene names 1 and 2 from blast file
        bName2 = lnLst[1]
        if oName1 == bName1 and oName2 == bName2:
            out.write(ln)
    oCnt += 1
    print oCnt



ortho.close()
blast.close()
out.close()
