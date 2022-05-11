#This script is designed to take an expression divergence file generated from
#cuffdiff and output a file that is a subset fo the input.  All lines that are
#output are lines where gene names are found in a reference input file.
#Created by David E. Hufnagel on Oct 25, 2012
import sys

inp = open(sys.argv[1])
ref = open(sys.argv[2])
out = open(sys.argv[3], "w")




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and take genes in the 1st col into a list
refLst = []
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        refLst.append(lineLst[0])

#Go through inp
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if lineLst[0] in refLst:
            out.write(line)



inp.close()
ref.close()
out.close()
