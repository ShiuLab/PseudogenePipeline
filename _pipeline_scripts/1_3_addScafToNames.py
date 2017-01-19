#This script is designed to add "scaffold_" to the chromosome number in column
#2 of a 4col file
#Designed by David E. Hufnagel on 5-25-2012

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")

for line in inp:
    lineLst = line.split("\t")
    newName = "scaffold_" + lineLst[1]
    newLine = "%s\t%s\t%s\t%s" % (lineLst[0], newName, lineLst[2], lineLst[3])
    out.write(newLine)



inp.close()
out.close()
