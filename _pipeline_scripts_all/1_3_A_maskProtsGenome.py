#This script is designed to mask a genome file's protein coding regions with
#"Ms" for "masked"
#Created by David E. Hufnagel on Aug 3, 2012
#WARNING: SITUATION SPECIFIC (Added an At to the fasta names)

import sys

fasta = open(sys.argv[1])    #input fasta format genome file
prot = open(sys.argv[2])     #input 4col proteins file
out = open(sys.argv[3], "w") #output masked fasta format genome file



def OrderCoords(coor1, coor2):
    if type(coor1) != int or type(coor2) != int:
        print "***TYPE ERROR***"
    
    if coor1 < coor2:
        small = str(coor1)
        big = str(coor2)
    elif coor2 < coor1:
        small = str(coor2)
        big = str(coor1)
    else:
        print "***ERROR HERE***"

    return small, big

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

def MaskRegion(seq, start, stop, newChar):  #start and stop should be indexes.  newChar is the character to mask with
    maskedChunk = (int(stop) - int(start) + 1) * newChar
    newSeq = seq[:start] + maskedChunk + seq[stop+1:]
    return newSeq
    



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#go through protein file and make dict with key: chromo val: [(startCoord, endCoord),]
protDict = {}
for line in prot:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[1], (lineLst[2],lineLst[3].strip()), protDict)

#go through fasta file and output masked sequences
currChromo = ""
for line in fasta:
    if not line.startswith("#"):
        if line.startswith(">"):
            #currChromo = "At" + line.strip().strip(">")  #for At specifically
            currChromo = line.strip().strip(">")
            out.write(line)
            cnt = 0
        else:
            #throw an error if there are multiple fasta sequence lines
            if cnt != 0:
                print "\n***error multiple seq lines per seq!***\n"
            cnt += 1

            #do the actual masking
            seq = line.strip()
            if currChromo in protDict:
                for tup in protDict[currChromo]:
                    startCoord, endCoord = OrderCoords(int(tup[0]), int(tup[1]))
                    startInd = int(startCoord) - 1
                    endInd = int(endCoord) - 1
                    seq = MaskRegion(seq, startInd, endInd, "M")

            #output new masked line
            out.write(seq + "\n")



fasta.close()
prot.close()
out.close()
