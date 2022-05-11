#This script is designed to extract pseudogenes from the MCScanX .tandem file
#containing pseudogenes and genes and output the info in a 4col format
#Created by David E. Hufnagel on May 9, 2013
import sys

four = open(sys.argv[1])     #input 4col pseudogenes file
tandem = open(sys.argv[2])   #input MCScanX .tandem file
out = open(sys.argv[3], "w") #output 4col tandem pseudogenes file



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through tandem and make a list of tandem pseudogenes
tandemLst = []
for line in tandem:
    if not line.startswith("#"):
        lineLst = line.strip().split(",")
        for item in lineLst:
            if item.startswith("Ps"):
                tandemLst.append(item)

#Go though four and output all lines containing tandem pseudogenes
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        psName = lineLst[0]
        if psName in tandemLst:
            out.write(line)
            


four.close()
tandem.close()
out.close()
