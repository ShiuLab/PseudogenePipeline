#This script is designed to compare a reciprical best match file and an
#orthologs file to see how much overlap there is
#Created by David E. Hufnagel on 3-9-2012

import sys

rbm = open(sys.argv[1])  #reciprocal best match file
orth = open(sys.argv[2]) #ortholog file



rbmDict = {}
for line in rbm:
    lineLst = line.split("\t")
    rbmDict[lineLst[0]] = lineLst[1][:-1]
    rbmDict[lineLst[1][:-1]] = lineLst[0]

goodCnt = 0 #overlapping matches
totalCnt = 0
for ln in orth:
    lnLst = ln.split("\t")
    x = lnLst[0]
    y = lnLst[1]
    if x in rbmDict:
        if rbmDict[x] == lnLst[1][:-1]:
            goodCnt += 1
    totalCnt += 1

bad = (totalCnt - goodCnt)
print
print "good: ", goodCnt, float(goodCnt)/totalCnt * 100 ,"%"
print "bad: ", bad, float(bad)/totalCnt * 100 ,"%"
print "total: ", totalCnt
print


rbm.close()
orth.close()
