#This script is designed to make a 2col homology file (for MCScanX) from a 4col
#pseudogene file with code names.
#Created by David E. Hufnagel on Mar 21, 2013
import sys

four = open(sys.argv[1])
out = open(sys.argv[2], "w")



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        pseuName = lineLst[0]
        geneName = "_".join(lineLst[0].split("_")[1:])
        newLine = "%s\t%s\n" % (pseuName, geneName)
        out.write(newLine)



four.close()
out.close()
