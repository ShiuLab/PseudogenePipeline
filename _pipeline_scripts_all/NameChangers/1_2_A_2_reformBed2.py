#This script was designed to take a .bed file and remove the decimals from the
#AT gene names.
#Created by David E. Hufnagel on 4-17-2012

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



def FixName(old):
    if old.startswith("AT"):
        new = old.split(".")[0]
        return new
    else:
        return old



for line in inp:
    lineLst = line.split("\t")
    newName1 = FixName(lineLst[0])
    newName2 = FixName(lineLst[1][:-1])
    newline = "%s\t%s\n" % (newName1, newName2)
    out.write(newline)



inp.close()
out.close()
