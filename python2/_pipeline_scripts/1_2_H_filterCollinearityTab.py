#This script is designed to take a full-length tabular (Guarav style) collinearity
#file and parse out certain species pairs.
#input format:
#G1     G2      BlockID|PairID|EVALUE
#Created by David E. Hufnagel on Jan 23, 2013


import sys

inp = open(sys.argv[1])      #input tabular collinearity file
pairs = sys.argv[2]          #pairs of species unique identifiers at the
                             #beginning of gene names.  pairs seperated by "_"
                             #and genes seperated by "," Ex: AL,Bra_AT,AT
out = open(sys.argv[3], "w") #output tabular collinearity file




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#go through pairs and put them into a list of tuples Ex: [(AT, Bra),(AT, AT)]
pairLst = []
for pair in pairs.split("_"):
    pairLst.append(tuple(pair.split(",")))
    
#go through inp and output lines that follow the criteria defined by pairs
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split()
        for pair in pairLst:
            if (lineLst[0].startswith(pair[0]) and lineLst[1].startswith(pair[1])) \
               or (lineLst[0].startswith(pair[1]) and lineLst[1].startswith(pair[0])):
                out.write(line)
    else:
        out.write(line)
        



inp.close()
out.close()
