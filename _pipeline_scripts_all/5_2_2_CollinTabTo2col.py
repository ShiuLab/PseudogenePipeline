#This script is designed to make a 2col pairs file with all gene pairs present
#in a Dave-style tabular .collinearity MCScanX file
#Created by David E. Hufnagel on April 10, 2013

import sys

inp = open(sys.argv[1])      #Input Dave-style tabular .collinearity MCScanX file
out = open(sys.argv[2], "w") #Output 2col pairs file



#Write the users command line prompt on the first line of the output file.
out.write("#python %s" % (" ".join(sys.argv)))

#Do the actual processing
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        for pair in lineLst[6:]:
            out.write("\n")
            out.write("\t".join(pair.split(";")))




inp.close()
out.close()
