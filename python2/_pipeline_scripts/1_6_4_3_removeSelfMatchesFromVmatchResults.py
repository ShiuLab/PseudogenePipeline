#This script is designed to remove self matches from Vmatch results files
#Created by David E. Hufnagel on July 17, 2012

import sys

inp = open(sys.argv[1])      #input Vmatch output file with self matches
out = open(sys.argv[2], "w") #output file without self matches





#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if not lineLst[1] == lineLst[5]:
            out.write(line)
            



inp.close()
out.close()
