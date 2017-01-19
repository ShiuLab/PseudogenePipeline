# This script is intended to make my .blast file names match my .bed file names.
# In this case I will simply remove the tenth decimal place from the
# A. Thaliana gene names.
# Created by David E. Hufnagel 10-25-2011

import sys

blast = open(sys.argv[1])    #the name of the .blast file to be processed
out = open(sys.argv[2], "w") #the name of the output

#print the unix command line that called this script
out.write('#python %s\n'%(' '.join(sys.argv)))

for line in blast:
    lineLst = line.split("\t")
    if not lineLst[0].startswith("AT"):
        print "ERROR HERE"
    name = lineLst[0].split(".")[0]
    newline = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (name, lineLst[1], lineLst[2], lineLst[3],\
                                lineLst[4], lineLst[5], lineLst[6], lineLst[7],\
                                lineLst[8], lineLst[9], lineLst[10], lineLst[11])
    out.write(newline)

blast.close()
out.close()
