#This script is designed to remove the "At" from the chromo names.  Ex: AtChr1
# --> Chr1 in a 4col file where the chromo is in the 2nd column.
#Created by David E. Hufnagel on Nov 27, 2012

import sys

inp = open(sys.argv[1])      #input 4col file with original chromo name
out = open(sys.argv[2], "w") #output 4col file with modified chromo name
spe = sys.argv[3]            #2 character species identifier at the beginning of the chromo


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        newName = lineLst[1].strip(spe)
        newLine = "%s\t%s\t%s\t%s\n" % (lineLst[0], newName, lineLst[2], lineLst[3])
        out.write(newLine)

inp.close()
out.close()
