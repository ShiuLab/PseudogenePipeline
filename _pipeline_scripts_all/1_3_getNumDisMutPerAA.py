#This script is designed to take a disable_count pseudogene file and outputing
#the number of disabling mutations / total AA in the format:
#name    numMut    totalAA    numMut/totalAA
#Created by David E. Hufnagel on June 4, 2013
import sys

inp = open(sys.argv[1])     #input .disable_count pseudogene file
ref = open(sys.argv[2])     #ref file with codeNames in the 1st col and bigName in the 2nd col
out = open(sys.argv[3],"w") #output .disMut file with information about disabling mutations



def OrderCoords(coor1, coor2):
    coor1 = int(coor1);coor2 = int(coor2)
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        print "***ORDERCOORDS ERROR HERE***"
        return coor1,coor2

    return small, big




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and make a dict of key: bigName val: codeName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        refDict[lineLst[1]] = lineLst[0]
        
#Go through inp, determine numMut, totalAA and ratio (numMut / totalAA) and output info
for line in inp:
    if not line.startswith("#python"):
        if line.startswith("#"):
            lineLst = line.strip().split(" ")
            small, big = OrderCoords(lineLst[2].split(":")[0].split("-")[0],\
                                lineLst[2].split(":")[0].split("-")[1])
            numAA = big - small + 1
            numDis = int(lineLst[4]) + int(lineLst[5]) + int(lineLst[6]) + int(lineLst[7])
            ratio = float(numDis) / float(numAA)
            codeName = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:]))[1:]
            psName = refDict[codeName]
            newLine = "%s\t%s\t%s\t%s\n" % (psName, numDis, numAA, ratio)
            out.write(newLine)
            




inp.close()
ref.close()
out.close()
