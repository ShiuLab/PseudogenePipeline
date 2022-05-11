#This script is designed to change the names in a fasta file from zfill(7) to
#zfill(6)
#Created by David E. Hufnagel on July 31, 2012
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])      #input fasta file with zfill(7) names
out = open(sys.argv[2], "w") #output fasta file with zfill(6) names



for line in inp:
    if not line.startswith("#"):
        if line.startswith(">"):
            num = line.split("Fake")[1]
            newNum = num[1:]
            newName = ">AtFake%s" % (newNum)
            out.write(newName)
        else:
            out.write(line)


inp.close()
out.close()
