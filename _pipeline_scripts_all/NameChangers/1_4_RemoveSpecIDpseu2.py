#This script is intended to remove species names from a fasta file
#Ex: ChlNC64A_1|scaffold_1 --> scaffold_1
#Created by David E. Hufnagel June 21, 2012

import sys

inp = open(sys.argv[1])      #The input fasta file
out = open(sys.argv[2], "w") #The output fasta file





for line in inp:
    if line.startswith(">"):
        newName = ">%s" % (line.split("|")[1])
        out.write(newName)
    else:
        out.write(line)





inp.close()
out.close()
