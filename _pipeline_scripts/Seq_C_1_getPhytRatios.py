#This script is designed to get the ratios of At-Br in at_br_block.csv
#Created by David E. Hufnagel on 4-10-2012
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

csv = open(sys.argv[1])     #at_br_block.csv

#go through csv
csvDict = {}
csv.readline()
for ln in csv:
    lnLst = ln.split(",")
    if lnLst[3] not in csvDict:
        csvDict[lnLst[3]] = [lnLst[4]]
    else:
        csvDict[lnLst[3]].append(lnLst[4])

print csvDict

#get ratios
oneCnt = 0
twoCnt = 0
threeCnt = 0
fourCnt = 0
for ln in csvDict.values():
    print ln
    if len(ln) == 1:
        oneCnt += 1
    elif len(ln) == 2:
        twoCnt += 1
    elif len(ln) == 3:
        threeCnt += 1
    elif len(ln) > 3:
        fourCnt += 1
    else:
        print "WHAT!!!!!"

total = float(oneCnt + twoCnt + threeCnt + fourCnt)
print "orthos At->Br"
print "1:1:  %s" % (oneCnt / total * 100)
print "1:2:  %s" % (twoCnt / total * 100)
print "1:3:  %s" % (threeCnt / total * 100)
print "1:4+: %s" % (fourCnt / total * 100)

            

csv.close()
