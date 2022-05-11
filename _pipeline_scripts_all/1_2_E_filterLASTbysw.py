#This script is designed to filter a .last output by sw score
#Designed by David E. Hufnagel on June 15, 2012

import sys

inp = open(sys.argv[1])      #input unfiltered .last file (.last)
out = open(sys.argv[2], "w") #output filtered .last file (.last.filt)
swThresh = int(sys.argv[3])  #Smith-Waterman score threshold





#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        sw = int(lineLst[0])
        if sw >= swThresh:
            out.write(line)
            #print line





inp.close()
out.close()
