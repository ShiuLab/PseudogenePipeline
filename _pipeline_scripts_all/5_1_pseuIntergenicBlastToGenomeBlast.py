#This script is designed to take a 4col file with intergenic genomic coordinates
#and make a 4col file with genomic coordinates.
#Created by David E. Hufnagel on Feb 12, 2013

import sys

blast = open(sys.argv[1])    #input 4col pseudogene file with intergenic coordinates
ref = open(sys.argv[2])      #input 4col intergenic coordinates file
out = open(sys.argv[3], "w") #output 4col pseudogene file with genomic coordinates




def OrderCoords(coor1, coor2):
    coor1 = int(coor1);coor2 = int(coor2)
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        print "***ERROR HERE***"
        return coor1,coor2

    return small, big



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#go through ref and make a dict of key: name val(chromo, startCoord, endCoord)
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        small, big = OrderCoords(int(lineLst[2]), int(lineLst[3]))
        refDict[lineLst[0]] = (lineLst[1], small, big)

#go through blast and output pseudogene info with new coordinates
for line in blast:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        small, big = OrderCoords(int(lineLst[8]), int(lineLst[9]))
        newStart = small + int(refDict[lineLst[1]][1]) - 1
        newEnd = big + int(refDict[lineLst[1]][1]) - 1
        chromo = refDict[lineLst[1]][0]
        mid = "\t".join(lineLst[2:8])#this middle chunk of the line from lineLst[2] to lineLst[7]        
        newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" % (lineLst[0], chromo, mid, newStart, newEnd, lineLst[10], lineLst[11])
        out.write(newLine)
        



blast.close()
ref.close()
out.close()
