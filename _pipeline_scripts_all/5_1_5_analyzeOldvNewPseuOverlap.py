#This script was designed to compare pseudogene files generated from the two
#pseudogene methods (with intergenic sequence (new) vs. with overlap pipeline and
#the original genome (old)).  This script takes an output file from Gaurav's
#overlap pipeline.
#Created by David E. Hufnagel on Jan 28, 2013
import sys

inp = open(sys.argv[1])      #input .onlyoverlap file
out = open(sys.argv[2], "w") #output info file



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

#gather info
sameCnt = 0; sameLst = []           #for cases where both things are the same
firstLongCnt = 0; firstLongLst = [] #for cases where the first thing is longer than the second on one end
secLongCnt = 0; secLongLst = []     #for cases where the second thing is longer than the first on one end
firstBigCnt = 0; firstBigLst = []   #for cases where the first thing is longer than the second on both ends
secBigCnt = 0; secBigLst = []       #for cases where the second thing is longer than the first on both 
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        smallFirst, bigFirst = OrderCoords(lineLst[2], lineLst[3])
        smallSec, bigSec = OrderCoords(lineLst[6], lineLst[7])
        lenFirst = bigFirst - smallFirst + 1
        lenSec = bigSec - smallSec + 1
        
        #for cases where both things are the same
        if lineLst[1] == lineLst[5] and ((lineLst[2] == lineLst[6] and lineLst[3] == lineLst[7]) or (lineLst[2] == lineLst[7] and lineLst[3] == lineLst[6])):
            sameLst.append(line)
            sameCnt += 1
        
        elif lineLst[1] == lineLst[5] and ((lineLst[2] == lineLst[6] or lineLst[2] == lineLst[7]) or (lineLst[3] == lineLst[6] or lineLst[3] == lineLst[7])):
            #for cases where the first thing is longer than the second on one end
            if lenFirst > lenSec:
                firstLongLst.append(line)
                firstLongCnt += 1
            #for cases where the second thing is longer than the first on one end
            elif lenSec > lenFirst:
                secLongLst.append(line)
                secLongCnt += 1

        else:
            #for cases where the first thing is longer than the second on both ends
            if lenFirst > lenSec:
                firstBigLst.append(line)
                firstBigCnt += 1
            #for cases where the second thing is longer than the first on both 
            elif lenSec > lenFirst:
                secBigLst.append(line)
                secBigCnt += 1

#output info
out.write("## same on both sides:                %s\n" % (sameCnt))
for line in sameLst:
    out.write(line)
out.write("## one side same, first side longer:  %s\n" % (firstLongCnt))
for line in firstLongLst:
    out.write(line)
out.write("## one side same, second side longer: %s\n" % (secLongCnt))
for line in secLongLst:
    out.write(line)
out.write("## first longer on both sides         %s\n" % (firstBigCnt))
for line in firstBigLst:
    out.write(line)
out.write("## second longer on both sides        %s\n" % (secBigCnt))
for line in secBigLst:
    out.write(line)
