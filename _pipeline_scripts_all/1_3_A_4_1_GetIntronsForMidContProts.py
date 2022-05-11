#This script is designed to get an intron dataset from a subset of the total
#number of proteings
#Created by David E. Hufnagel on Aug 16, 2012

import sys

introns = open(sys.argv[1])  #input introns 4col file
midProts = open(sys.argv[2]) #input proteins to be used 4col file
out = open(sys.argv[3], "w") #output intron 4col file with only introns from proteins in midProts




def ExtractColFromFileIntoList(ind, fd):  #fd = file directory
    itemLst = []
    for line in fd:
        if not line.startswith("#"):
            lineLst = line.split("\t")
            itemLst.append(lineLst[ind])

    return itemLst




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#go through midProts and make a list of names
protLst = ExtractColFromFileIntoList(0, midProts)

#go through introns and output lines where the protein names are in protLst
for line in introns:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        prot = lineLst[0].split("Intron")[0]
        if prot in protLst:
            out.write(line)





introns.close()
midProts.close()
out.close()
