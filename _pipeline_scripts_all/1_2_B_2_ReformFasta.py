#This script is intended to add the species ID to the names of a fasta file
#WARNING: Highy situation specific
#Created by David E. Hufnagel on 2-17-2012

import sys

inp = open(sys.argv[1])       #The input fasta file
out = open(sys.argv[2], "w")  #The new fasta file

for line in inp:
    if line.startswith(">"): #the name lines
        newline = line
        if line[1:3] == "Al":
            pass
        elif line[1:3] == "AT":
            newline = ">At|" + line[1:]
        elif line[1] == "B":
            newline = ">Br|" + line[1:]
        elif line[1] == "R":
            newline = ">Rr|" + line[1:]
        else:
            print "***ERROR HERE!***"
            #print line
        print newline
        out.write(newline)
    else:
        out.write(line)
    
    



inp.close()
out.close()
