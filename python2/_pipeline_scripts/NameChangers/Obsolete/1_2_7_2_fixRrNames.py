#This script is designed to change names in the second column from
#RrC7747_p1_1.000 --> RrC7747_p1
#Created by David E. Hufnagel on 3-16-2012

import sys

inp = open(sys.argv[1])      #input file
out = open(sys.argv[2], "w") #output file

for line in inp:
    lineLst = line.split("\t")
    newName = "_".join(lineLst[1][:-1].split("_")[:-1])
    newLine = "%s\t%s\n" % (lineLst[0], newName)
    out.write(newLine)

inp.close()
out.close()
