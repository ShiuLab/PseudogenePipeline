#This script is designed to parse an edgeR output file by maximum FDR.
#Created by David E. Hufnagel on Nov

import sys

inp = open(sys.argv[1])          #input FDR file
out = open(sys.argv[2], "w")     #output parsed FDR file
FDRthresh = float(sys.argv[3])   #FDR threshold
#logFC threshold is assumed to be smaller or equal to -1 or greater or equal to 1


def RemoveBlancs(theList):
    newLst = []
    for item in theList:
        if item != "":
            newLst.append(item)

    return newLst



out.write(inp.readline())
for line in inp:
    lineLst = line.strip("\n").split(" ")
    lineLst = RemoveBlancs(lineLst)

    if float(lineLst[-1]) <= FDRthresh and (float(lineLst[1]) <= -1.0 or float(lineLst[1]) >= 1.0):
        out.write(line)




inp.close()
out.close()
