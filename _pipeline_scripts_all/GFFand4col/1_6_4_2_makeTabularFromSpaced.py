#This script is designed to replace all spaces in a tabular file with one \t
#Created by David E. Hufnagel on July 17, 2012

import sys

inp = open(sys.argv[1])      #input spaced out file
out = open(sys.argv[2], "w") #output modified tabular file





for line in inp:
    if not line.startswith("#"):
        lineLst = line.split(" ")
        newLineLst = []
        for item in lineLst:
            if item != "":
                newLineLst.append(item)
        newLine = "\t".join(newLineLst)
        out.write(newLine)
    else:
        out.write(line)
            




inp.close()
out.close()
