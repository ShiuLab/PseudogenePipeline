#This script is part of a series designed to test the rates of putative genomic
#missasembly. More specifically, it takes the distance between the EST match and
#the end of the EST (Y), the distance between the matching region on the genome
#and the end of the contig (X) and the average intron length for the species of
#interest (Z) to calculate whether Y <= X + Z  If true, add a missasemby
#evidence to the contig, with a bit more complexity.  Part III is designed to
#take the part2 ests output file and condense/convert it to a contig based, final
##output file in the following format:
#contigName  estNum  misNum  goodNum  ambNum  estLst
#with one line per est
#Created by David E. Hufnagel on August 21, 2012

import sys

inp = open(sys.argv[1])      #input part 2 est match file
out = open(sys.argv[2], "w") #output part 3 final contig est matches file




def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)



        
#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
out.write("#contigName\testNum\tmisNum\tgoodNum\tambNum\testLst\n")

#Get input into dict of key: contigName val: [(estName, isMisassembled?),]
inpDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[1], (lineLst[0], lineLst[2].strip()), inpDict)

#Go through dict, process and output info
for contigName in inpDict:
    #Calculate total misassembly, total good, total ambiguous and total ests mapped
    estCnt = len(inpDict[contigName])
    misCnt = 0
    goodCnt = 0
    ambCnt = 0
    estLst = []
    for est in inpDict[contigName]:
        estLst.append(est[0])
        if est[1] == "misassembled":
            misCnt += 1
        elif est[1] == "good":
            goodCnt += 1
        elif est[1] == "ambiguous":
            ambCnt += 1
        else:
            print "ERROR HERE!"

    #Output info
    estStr = ",".join(estLst)
    newLine = "%s\t%s\t%s\t%s\t%s\t%s\n" % (contigName, estCnt, misCnt, goodCnt, ambCnt, estStr)
    out.write(newLine)
    
     


inp.close()
out.close()
