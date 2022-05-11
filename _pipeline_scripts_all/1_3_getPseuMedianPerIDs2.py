#This script is designed to acquire %IDs for all pseudoexons, calculate the
#median and output the information in a file with the format:
#chromo;start-stop  exonCoord1;exonCoord2...  %ID1;%ID2...  median%ID
#exonCoord format:
#coord1-coord2
#Created by David E. Hufnagel on May 8, 2013
#WARNING DOESN'T WORK!!! FOR SOME REASON THERE IS AN OUTPUT LINE FOR EACH
#PSEUDOEXON INSTEAD OF EACH PSEUDOGENE
import sys

exon = open(sys.argv[1])     #input merged pseudoexon file (.PE_I457.PS1)
blast = open(sys.argv[2])    #input blast file used to find the pseudogenes
out = open(sys.argv[3], "w") #output file


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

#Take all pseudoexons and get the coordinates for the pseudogene in the form:
#startCoord;stopCoord
def GetStartStop(group):
    #extract all numbers into one list
    numLst = []
    for pair in eval(group):
        for num in pair:
            numLst.append(num)

    #take min and max and output the result
    mini = min(numLst)
    maxi = max(numLst)
    startStop = "%s-%s" % (str(mini), str(maxi))

    return startStop

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

#gets the median of a list
def median(listx):
    #Make all items floats and sort the list
    newLst = []
    for item in listx:
        newLst.append(float(item))
    newLst.sort()

    #Determine if the length is an odd number
    if len(newLst) % 2:
        isOdd = True
    else:
        isOdd = False

    #Find the actual median
    if isOdd == True:
        median = float(newLst[len(newLst)/2])
    else:
        leftMid = newLst[(len(newLst)/2)-1]
        rightMid = newLst[len(newLst)/2]
        median = float(leftMid + rightMid) / 2.0

    return median
    
#Go through exon and make a dict of key: chromo;gene;ExStart-ExStop val: chromo;gene;PsStart-PsStop
exonDict = {}
for line in exon:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        chromo = lineLst[0]
        gene = lineLst[1]
        startStop = GetStartStop(lineLst[2])
        psName = "%s;%s;%s" % (chromo, gene, startStop)
        for ex in eval(lineLst[2]):
            key = "%s;%s;%s-%s" % (chromo, gene, ex[0], ex[1])
            exonDict[key] = psName

#Go through blast and make a dict of key: psName val: [%ID1,%ID2...]
idDict = {}
for line in blast:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        chromo = lineLst[1]
        gene = lineLst[0]
        subjStart, subjStop = OrderCoords(lineLst[8], lineLst[9])
        perID = lineLst[2]
        exName = "%s;%s;%s-%s" % (chromo, gene, subjStart, subjStop)
        try:
            psName = exonDict[exName]
        except KeyError:
            continue
        SaveIntoDict(psName, perID, idDict)

        if chromo == "AthaChr1" and gene == "AT1G01010":
            print lineLst
            print exName
            print psName
            print idDict
            print 
        
#Go through blastDict, calculate median %ID and output info
for psName in idDict:
    ids = idDict[psName]
    medID = median(ids)
    out.write("%s\n" % (medID))



exon.close()
blast.close()
out.close()
