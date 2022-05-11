#This script is designed to Remove lines from a fasta file that are in a bad
#names input 4col file.
#Created by David E. Hufnagel on July 27, 2012

import sys

inp = open(sys.argv[1])      #input fasta file
bad = open(sys.argv[2])      #input 4col file with the bad names
out = open(sys.argv[3], "w") #output fasta file without bad name lines





#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#input bad into a list of names
badLst = []
for line in bad:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        name = lineLst[1]
        badLst.append(name)

#go through inp and keep only the lines without names in badLst
doKeep = True
for line in inp:
    if not line.startswith("#"):
        if line.startswith(">"):
            name = line.strip(">").strip().split(" ")[0]
            if name in badLst:
                doKeep = False
            else:
                doKeep = True
                out.write(line)
        else:
            if doKeep == True:
                out.write(line)
        





inp.close()
bad.close()
out.close()
