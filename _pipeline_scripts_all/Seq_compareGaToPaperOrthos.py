#This script is designed to compare AtBrOrthosBrPaper.stuff to
#AtBrOrthosGaurav2.xinp to see how many pairs are the same
#Designed by David E. Hufnagel on 4-19-2012
#WARNING HIGHLY SITUATION SPECIFIC

import sys

wang = open(sys.argv[1])     #AtBrOrthosBrPaper.stuff
moghe = open(sys.argv[2])    #AtBrOrthosGaurav2.xinp
#out = open(sys.argv[3], "w") #output file


def SaveIntoDict(gene1, gene2, dictX):
    if gene1 not in dictX:
        dictX[gene1] = [gene2]
    else:
        dictX[gene1].append(gene2)
    


#writes the users command line prompt on the first line of the output file.
out.write('#python %s\n'%(' '.join(sys.argv)))

#Go through wang and make a dict of At-Br pairs.  (and count total pairs)
wangDict = {}
wangCnt = 0
for line in wang:
    lineLst = line.split("\t")
    if lineLst[1] != "-" and \
       (lineLst[4] != "-" or lineLst[5] != "-" or lineLst[6][:-1] != "-"):
        wangCnt += 1
        if (lineLst[4] != "-"):
            SaveIntoDict(lineLst[1], lineLst[4], wangDict)
        if (lineLst[5] != "-"):
            SaveIntoDict(lineLst[1], lineLst[5], wangDict)
        if (lineLst[6][:-1] != "-"):
            SaveIntoDict(lineLst[1], lineLst[6][:-1], wangDict)

#Go through moghe and make a dict of At-Br pairs. (and count total pairs)
mogheDict = {}
mogheCnt = 0
for ln in moghe:
    mogheCnt += 1
    lnLst = ln.split("\t")
    if lnLst[0] not in mogheDict:
        mogheDict[lnLst[0]] = [lnLst[1][:-1]]
    else:
        mogheDict[lnLst[0]].append(lnLst[1][:-1])

#Compare wang and moghe dicts and count pair overlap.

intCnt = 0
for pair in wangDict:
    if pair in mogheDict:
        for x in wangDict[pair]:
            for y in mogheDict[pair]:
                if x == y:
                    intCnt += 1

print wangCnt
print mogheCnt
print intCnt




wang.close()
moghe.close()
