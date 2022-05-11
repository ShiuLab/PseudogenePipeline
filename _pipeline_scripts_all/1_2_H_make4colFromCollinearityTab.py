#This script was designed to take a tabular (Gaurav style) colinearity file from
#MCSCanX and make a 4COL FILE OF SYNTENIC BLOCKS from it.
#Created by David E. Hufnagel on Jan 23, 2013

import sys

blocks = open(sys.argv[1])    #input tabular (Gaurav style) colinearity file
four = open(sys.argv[2])      #input 4col genes file (must contain all genes in blocks)
out1 = open(sys.argv[3], "w") #1st output 4col syntenic blocks file
out2 = open(sys.argv[4], "w") #2nd output 4col syntenic blocks file



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



#Write the users command line prompt on the first line of the output files.
out1.write("#python %s\n" % (" ".join(sys.argv)))
out2.write("#python %s\n" % (" ".join(sys.argv)))

#Go through 4col file and make refDict of key = geneName val = (chromo, startCoord, endCoord)
refDict = {}
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        small, big = OrderCoords(lineLst[2], lineLst[3])
        refDict[lineLst[0]] = (lineLst[1], small, big)

#Go through collinearity file, get min and max coordinates for all genes in blocks and ouput info       
start1 = 999999999 #the start coordinate for the syntenic block in the 1st chromosome
stop1 = 0          #the end coordinate for the syntenic block in the 1st chromosome
start2 = 999999999 #the start coordinate for the syntenic block in the 2nd chromosome
stop2 = 0          #the end coordinate for the syntenic block in the 2nd chromosome
lastBlockID = ""
chromo1 = ""
chromo2 = ""
blockName = ""
cnt = 0
for line in blocks:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        blockID = lineLst[2].split("|")[0]
        if not blockID == lastBlockID:
            chromo1 = refDict[lineLst[0]][0]
            chromo2 = refDict[lineLst[1]][0]
            if cnt != 0:
                blockName = "block_%s_%s|%s" % (lastBlockID, lastChromo1, lastChromo2)
                newLine1 = "%s_%s\t%s\t%s\t%s\n" % (lastChromo1, blockName, lastChromo1, start1, stop1)
                newLine2 = "%s_%s\t%s\t%s\t%s\n" % (lastChromo2, blockName, lastChromo2, start2, stop2)
                out1.write(newLine1)
                out2.write(newLine2)
                start1 = 999999999
                stop1 = 0
                start2 = 999999999
                stop2 = 0
            cnt += 1
        tempStart1 =  int(refDict[lineLst[0]][1])
        if tempStart1 < start1:
            start1 = tempStart1
        tempStart2 = int(refDict[lineLst[1]][1])
        if tempStart2 < start2:
            start2 = tempStart2
        tempStop1 = int(refDict[lineLst[0]][2])
        if tempStop1 > stop1:
            stop1 = tempStop1
        tempStop2 = int(refDict[lineLst[1]][2])
        if tempStop2 > stop2:
            stop2 = tempStop2

        lastBlockID = blockID
        lastChromo1 = chromo1
        lastChromo2 = chromo2
else:
    blockName = "block_%s_%s|%s" % (lastBlockID, lastChromo1, lastChromo2)
    newLine1 = "%s_%s\t%s\t%s\t%s\n" % (lastChromo1, blockName, lastChromo1, start1, stop1)
    newLine2 = "%s_%s\t%s\t%s\t%s\n" % (lastChromo2, blockName, lastChromo2, start2, stop2)
    out1.write(newLine1)
    out2.write(newLine2)
                



blocks.close()
four.close()
out1.close()
out2.close()

