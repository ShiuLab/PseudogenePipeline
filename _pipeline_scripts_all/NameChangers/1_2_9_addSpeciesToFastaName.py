#This script was intended to add a species identifier to every name in a
#fasta file
#created by David E. Hufnagel on 2-10-2012

import sys

inp = open(sys.argv[1])  #the name of the fasta input file
spe = sys.argv[2]        #the species identifier
out = open(sys.argv[3], "w") #the name of the output file

for line in inp:
    if line.startswith(">"):  #on title lines
        newLine = ">" + spe + "|" + line[1:]
        out.write(newLine)
    else:
        out.write(line)

inp.close()
out.close()
