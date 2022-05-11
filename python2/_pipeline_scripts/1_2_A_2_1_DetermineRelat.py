#This script is designed to take a MCScanX collinearity file and determine the
#kind of relationships between the genes in the blocks (1:1 vs. 1:2 vs 1:3 etc)
#Created by David E. Hufnagel on 3-21-2012

import sys

col = open(sys.argv[1])      #input collinearity file
out = open(sys.argv[2], "w") #output file to contain the data



def SkipHash(fd):
    for L in fd:
        if not L.startswith("#"):
            break



#write the users command line prompt on the first line of the output file.
out.write('#python %s\n'%(' '.join(sys.argv)))
#skip past the title in the collinearity file
SkipHash(col)

#go through the collinearity file and make the dict and the inverse dict
bigDict = {}
invDict = {}
for line in col:
    if not line.startswith("##"):
        lineLst = line.split("\t")
        #for bigDict
        if lineLst[1] not in bigDict:
            bigDict[lineLst[1]] = [lineLst[2]]
        else:
            bigDict[lineLst[1]].append(lineLst[2])

        #for invDict
        if lineLst[2] not in invDict:
            invDict[lineLst[2]] = [lineLst[1]]
        else:
            invDict[lineLst[2]].append(lineLst[1])

##        if lineLst[1] == "AT5G06600":
##            print lineLst
##            print bigDict[lineLst[1]]

#go through dictionaries to count relationships
#forward
oneCnt1 = 0
twoCnt1 = 0
threeCnt1 = 0
moreCnt1 = 0
for key in bigDict:
    if len(bigDict[key]) == 1:
        oneCnt1 += 1
    elif len(bigDict[key]) == 2:
        twoCnt1 += 1
    elif len(bigDict[key]) == 3:
        threeCnt1 += 1
    elif len(bigDict[key]) >= 4:
        moreCnt1 += 1
        print bigDict[key]
    else:
        print "\n\nERROR1!!!\n\n"

#inverse
oneCnt2 = 0
twoCnt2 = 0
threeCnt2 = 0
moreCnt2 = 0
for k in invDict:
    if len(invDict[k]) == 1:
        oneCnt2 += 1
    elif len(invDict[k]) == 2:
        twoCnt2 += 1
    elif len(invDict[k]) == 3:
        threeCnt2 += 1
    elif len(invDict[k]) >= 4:
        moreCnt2 += 1
    else:
        print "\n\nERROR1!!!\n\n"

#output data
out.write("1:1 col1:  %s\n" % (oneCnt1))
out.write("1:2 col1:  %s\n" % (twoCnt1))
out.write("1:3 col1:  %s\n" % (threeCnt1))
out.write("1:4+ col1: %s\n" % (moreCnt1))
out.write("\n")
out.write("1:1 col2:  %s\n" % (oneCnt2))
out.write("1:2 col2:  %s\n" % (twoCnt2))
out.write("1:3 col2:  %s\n" % (threeCnt2))
out.write("1:4+ col2: %s\n" % (moreCnt2))




col.close()
out.close()
