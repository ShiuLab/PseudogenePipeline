#This script is designed to remove all funky stuff from pseudogene amino acid
#files by simply deleting them.  Defenition of funky stuff: *, -, | and /
#Created by David E. Hufnagel on Dec 19, 2012

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        if line.startswith(">"):
            out.write(line)
        else:
            outSeq = ""
            badLst = "*-|/"
            for char in line.strip():
                if not char in badLst:
                    outSeq += char

            out.write(outSeq)
            out.write("\n")
                



inp.close()
out.close()
