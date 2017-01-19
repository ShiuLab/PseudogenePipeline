#This script is designed to get information around the overlapping of genome
#alignments from the LAST output.
#Created by David E. Hufnagel on June 14, 2012



import sys

inp = open(sys.argv[1])      #LAST output for processing
#lenT= int(sys.argv[2])       #length threshold
swT = int(sys.argv[2])       #Smith-Waterman score threshold
out = open(sys.argv[3], "w") #output file with information





def SaveIntoDict(gene1, gene2, dictX):
    if gene1 not in dictX:
        dictX[gene1] = [gene2]
    else:
        dictX[gene1].append(gene2)

#Not very useful, but kept arbitrarily
def GetRatio (val):
    chro = key.split(":")[0]
    ratio = float(val) / float(chroDict[chro]) * 100
    return chro, ratio




        
#go through input and make dictionary of alignments
bigDict = {}
chroDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        swScore = int(lineLst[0])
        
        if swScore >= swT:
            refName = lineLst[1]
            mainName = lineLst[6]
            mainLen = lineLst[8]
            mainChLen = lineLst[10]
            SaveIntoDict(mainName, (refName, mainLen), bigDict)
            if mainName not in chroDict:
                chroDict[mainName] = mainChLen

#go through bigDict and condense it so that the keys are refChromo specific.
condDict = {}
for mainChromo in bigDict:
    tuples = bigDict[mainChromo]
    for tup in tuples:
        newKey = str(mainChromo) + ":" + tup[0]
        SaveIntoDict(newKey, int(tup[1]), condDict)

#go throught condDict and condense the value list into one value by taking
#the sum of the alignment lengths.
condDict2 = {}
for key in condDict:
    condDict2[key] = sum(condDict[key])

#output condDict2 info and make a dictioanary for the total alignment values    
totalDict = {}  #holds the total alignment for a main chromosome
for key in condDict2:
    val = condDict2[key]
    #get % alignments
    chro, ratio = GetRatio(val)
    SaveIntoDict(chro, ratio, totalDict)

#go through totalDict and send total alignment values into totalDict2
totalDict2 = {}
for key in totalDict:
    totalDict2[key] = sum(totalDict[key])

#go through totalDict2 and output everything
for totKey in totalDict2:
    totVal = totalDict2[totKey]

    #output cumulative valueus (Ex: scaffold_2)
    newLine = "%s\t%s\n" % (totKey, totVal)
    out.write(newLine)

    #output specific values (Ex: scaffold_2:Chr1)
    for key in condDict2:
        if key.split(":")[0] == totKey:
            val = condDict2[key]
            chro, ratio = GetRatio(val)
            newLine = "%s\t%s\n" % (key, ratio)
            out.write(newLine)



inp.close()
out.close()
