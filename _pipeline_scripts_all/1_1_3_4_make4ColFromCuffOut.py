#This script is designed to make a 4 col file from the Cufflinks output
#genes.fpkm_tracking file in the format:
#GeneName  GeneCoords  FPKM  FPKMConfidenceInterval
#Created by David E. Hufnagel on Oct 18, 2012

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")




cnt = 0
for line in inp:
    #on the first input line make an output info line
    
    if cnt == 0:
        newLine = "#GeneName\tGeneCoords\tFPKM\tFPKMConfidenceInterval\n"
        out.write(newLine)
    else:
        lineLst = line.split("\t")
        newLine = "%s\t%s\t%s\t%s|%s\n" % (lineLst[0], lineLst[6], lineLst[9], lineLst[10], lineLst[11])
        out.write(newLine)
    cnt += 1
