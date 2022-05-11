#This script is designed to make a .size file from a 4 col file
#Created by David E. Hufnagel on June 26, 2012

import sys

inp = open(sys.argv[1])      #the input 4col file
out = open(sys.argv[2], "w") #the output .size file





#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        size = abs(int(lineLst[3]) - int(lineLst[2])) + 1
        newLine = "%s\t%d\n" % (lineLst[0], size)
        out.write(newLine)




inp.close()
out.close()
