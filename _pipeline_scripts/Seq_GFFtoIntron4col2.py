#This script is designed to make an intron 4col file from a gff file with start
#and stop codons and exons. (Version 1 has genes, UTRs and CDSs).  Also
#this version determines the longest transcript itself.
#In short this expects a .gff format whereas the first version expects a .gff3 format
#Created by David E. Hufnagel on May 2, 2013
import sys, os

gff = open(sys.argv[1])      #input gff file
out = open(sys.argv[2], "w") #output intron 4col file
spe = sys.argv[3]            #species identifier Ex: Atha



def OrderCoords(coor1, coor2):
    coor1 = int(coor1);coor2 = int(coor2)
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        print "***MINOR ERROR HERE***"
        return coor1,coor2

    return small, big

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

def SaveIntoCntDict(key, val, dictX):
    val = int(val)
    if key not in dictX:
        dictX[key] = val
    else:
        dictX[key] += val



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through the gff and make a dict of key: geneName val: [(coordLen, transcriptName;exonNum),] #also exonName=transcriptName;exonNum
#and also a dict for the start and stop codons of key: geneName val: geneStart/geneStop (respectively)
#and a dict for the chromo key: geneName val: chromo
exonLenDict = {}
startDict = {}
stopDict = {}
chromoDict = {}
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        geneName = lineLst[8].split("name")[1].split(";")[0].strip()[1:-1].replace(".","_")

        chromo = lineLst[0]
        chromoDict[geneName] = chromo
        
        small, big = OrderCoords(lineLst[3], lineLst[4])
        coordLen = big - small + 1
        if lineLst[2] == "exon":
            transcriptName = geneName + "_trID" + lineLst[8].split("transcriptId")[1].split(";")[0].strip()#[1:-1]
            #save info into dict without exon number info
            SaveIntoDict(geneName, (coordLen, transcriptName), exonLenDict)
            #add exon number info
            exonLenDict[geneName]
            #sys.exit()
        elif lineLst[2] == "start_codon":
            startDict[geneName] = small
        elif lineLst[2] == "stop_codon":
            stopDict[geneName] = big

#Determine longest transcript by adding up CDSs and put the info into a list of transcript names
longestTranscripts = []
for geneName in exonLenDict:
    #add up CDSs
    CDSs = exonLenDict[geneName]
    CDSdict = {} #a temporary dict with key: transcriptName val: currentCountOfTranscriptLen
    for CDS in CDSs:
        SaveIntoCntDict(CDS[1], CDS[0], CDSdict)

    #Determine the longest transcript
    currMax = ("x", -1) #the current longest transcript in form (name, length)
    for transc in CDSdict:
        if CDSdict[transc] > currMax[1]:
            currMax = (transc, CDSdict[transc])

    #Put the transcript name into the list
    longestTranscripts.append(currMax[0])

#pull exons from the longest transcript out of the gff and make exon2LenDict with key: geneName val: (start, stop transcriptName;exonNum)
exon2LenDict = {}
gff.seek(0)
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        geneName = lineLst[8].split("name")[1].split(";")[0].strip()[1:-1].replace(".","_")
        small, big = OrderCoords(lineLst[3], lineLst[4])
        coordLen = big - small + 1
        if lineLst[2] == "exon":
            transcriptName = geneName + "_trID" + lineLst[8].split("transcriptId")[1].split(";")[0].strip()#[1:-1]

            #if transcriptName in longestTranscripts output info into exon2LenDict
            if transcriptName in longestTranscripts:
                SaveIntoDict(geneName, (small, big, transcriptName), exon2LenDict)

#Go through startDict and use the other dicts to determine intron lengths and
#Output the info using code names
for geneName in startDict:
    #gather basic info
    try:
        start = int(startDict[geneName]) #if missing a start codon, skip the gene
        stop = int(stopDict[geneName])   #if missing a stop codon, skip the gene
    except KeyError:
        continue
    
    chromo = chromoDict[geneName]
    smallSS, bigSS = OrderCoords(start, stop)
    exons = exon2LenDict[geneName]

    #make a list of tuples of exon coords and sort it by start coords
    exonLst = []
    for exon in exons:
        small, big = OrderCoords(exon[0], exon[1])
        exonLst.append((small, big))
    exonLst.sort()

    #remove overlapping features and make a merged exon list
    lastExon = (0, -1)
    exonLst2 = []  #the merged exon list
    for exon in exonLst:
        #overlapping features |---|
        #                       |---|
        if lastExon[1] >= exon[0] and lastExon[1] < exon[1]:
            lastExon = (lastExon[0],exon[1]) 

        #enveloping features |---------|
        #                       |---|
        elif lastExon[1] > exon[1]:
            continue #do nothing, but skip this exon

        #seperate features |---| |-----|
        else:
            #output lastExon info
            if not lastExon == (0, -1):
                #the first exon
                exonLst2.append(lastExon)
            lastExon = exon
    else:
        exonLst2.append(lastExon)

    #extract all tuples from the exonLst into a simple list of ordered coords
    exonLst3 = []
    for exon in exonLst2:
        exonLst3.append(exon[0])
        exonLst3.append(exon[1])

    #remove the first and last items from the list and make new tuples from what's left
    intronLst = []
    cnt = 1
    for coord in exonLst3[1:-1]:
        if not cnt % 2:
            intronLst[-1] = (intronLst[-1], coord-1)
        else:
            intronLst.append(coord+1)
        cnt += 1

    #output info
    codeNum = 1
    for intron in intronLst:
        codeName = "%s__%s__Intron__%s" % (spe, geneName, codeNum)
        newLine = "%s\t%s\t%s\t%s\n" % (codeName, chromo, intron[0], intron[1])
        codeNum += 1
        out.write(newLine)



gff.close()
out.close()





