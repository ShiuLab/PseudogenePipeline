#This script is designed to acquire %IDs for all pseudoexons, calculate the
#median and output the information in a file with the format:
#psName  exonCoord1;exonCoord2... %ID1;%ID2... median%ID
#exonCoord format:
#coord1-coord2
#Created by David E. Hufnagel on May 8, 2013
#WARNING: DOES NOT WORK.  DISCOVERED SOME DIFFICULTIES IN CONNECTING PSEUDOGENES
#WITH THEIR PSEUDOEXONS.  WILL DO A SIMPLER ALGORITHM WITHOUT NAMES IN VERSION 2
import sys

pseu = open(sys.argv[1])     #input 4col filtered pseudogens file
exon = open(sys.argv[2])     #input merged pseudoexon file (.PE_I457.PS1)
blast = open(sys.argv[3])    #input blast file used to find the pseudogenes
out = open(sys.argv[4], "w") #output file


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
    startStop = "%s;%s" % (str(mini), str(maxi))

    return startStop

#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through pseu and make a dict of key: chromo;start;stop)  val: pseuName
pseuDict = {}
for line in pseu:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        pseuName = lineLst[0]
        chromo = lineLst[1]
        start, stop = OrderCoords(lineLst[2], lineLst[3])
        pseuDict["%s;%s;%s" % (chromo, start, stop)] = pseuName
print pseuDict
        #sys.exit()

print len(pseuDict)

#Go through exon and make a dict of key: pseuName val: (chromo, [[start, stop],...])
exonDict = {}
for line in exon:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        chromo = lineLst[0]
        startStop = GetStartStop(lineLst[2])
        try:
            pseuName = pseuDict["%s;%s" % (chromo, startStop)]
        except KeyError:
           # print lineLst
            continue

        exonDict[pseuName] = (chromo, eval(lineLst[2]))
        #print lineLst
        #print chromo
        #print startStop
        #print pseuName
        #print

print len(exonDict)

#Go through blast and make a dict of key: pseuName val: [%ID1,%ID2...]

#Go through blastDict, calculate median %ID and output info




pseu.close()
exon.close()
blast.close()
out.close()
