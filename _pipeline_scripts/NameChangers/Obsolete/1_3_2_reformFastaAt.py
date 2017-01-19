#this script is designed to remove the decimals from At names in a fasta file
#Created by: David E. Hufnagel on 5-7-2012

import sys

inp = open(sys.argv[1])      #input file
out = open(sys.argv[2], "w") #output file

#print the unix command line that called this script
out.write('#python %s\n'%(' '.join(sys.argv)))



def FixName(old):
    if old.startswith("AT"):
        new = old.split(".")[0]
        return new
    else:
        return old



for line in inp:
    if line.startswith(">"):
        newName = FixName(line[1:])
        newLine = ">%s\n" % (newName)
        out.write(newLine)
    else:
        out.write(line)




inp.close()
out.close()
