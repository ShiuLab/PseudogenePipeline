#This script was designed to take one input file with the Ks rates between
#retained duplicates and calculate the timing of duplication based on the
#formula T= Ks / 2*.oo7
#Created by David E. Hufnagel on Dec 13, 2012

import sys

inp = open(sys.argv[1])      #input file with Ka/Ks analysis info
out = open(sys.argv[2], "w") #output file with col1: name1;name2col2: T




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        names = "%s;%s" % (lineLst[0], lineLst[1])
        T = float(lineLst[5]) / (2.0*.007)
        newLine = "%s\t%s\n" % (names, T)
        out.write(newLine)




inp.close()
out.close()
