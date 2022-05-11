#This script was designed to remove the decimals that follow the gene names
#Designed by David E. Hufnagel on 3-9-2012

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")

for line in inp:
    lineLst = line.split("\t")
    old = lineLst[1][:-1]
    temp = old.split("_")
    new = "_".join(temp[:-1])
    newLine = "%s\t%s\n" % (lineLst[0], new)
    out.write(newLine)




inp.close()
out.close()
