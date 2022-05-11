#This script is designed for changing the names within a fasta file
#Created by: David E. Hufnagel

import sys

FD = sys.argv[1]
OUT = sys.argv[2]
fd = open(FD)
out = open(OUT, "w")

#print the unix command line that called this script
print out.write('#python %s\n'%(' '.join(sys.argv)))

for line in fd:
    if line[0] == ">":
    #take out all the junk
    #WARNING POOR CODE
        name = line[1:].strip()
        nameLst = name.split("|")
        name = nameLst[-1]
        
        nameLst = name.split(",")
        name = nameLst[0]
        
        nameLst = name.split(" ")
        name = nameLst[-1]
        
    #put back the names
        out.write(">" + name + "\n")
    else:
        out.write(line)
        
fd.close()
out.close()
print "Done!"


