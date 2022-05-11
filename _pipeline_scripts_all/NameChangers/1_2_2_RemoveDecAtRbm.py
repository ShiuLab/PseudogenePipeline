# This script was intended to remove the decimals from A. Thaliana gene names
# in rbm files.  Because I didn't need lines 3 and 4 I left them out of the new
# rbm file.
#Created by David E. Hufnagel 11-04-2011

import sys

rbm = open(sys.argv[1])      #rbm file for input
out = open(sys.argv[2], "w") #output file

#print the unix command line that called this script
out.write('#python %s\n'%(' '.join(sys.argv)))

for line in rbm:
    lineLst = line.split("\t")
    if not lineLst[0].startswith("AT") or lineLst[0].startswith("#"):
        print "ERROR HERE"
    name = lineLst[0].split(".")[0]
    newline = "%s\t%s\n" % (name, lineLst[1])
    out.write(newline)

rbm.close()
out.close()
