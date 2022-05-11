#This script is designed to add pseudogenes from a pseudogene GFF file to genes
#in an ordered gene GFF file.
#Created by David E. Hufnagel on July 16, 2013
#WARNING: SITUATION SPECIFIC
import sys

gene = open(sys.argv[1])     #input gene gff file
pseu = open(sys.argv[2])     #input pseudogene gff file
out = open(sys.argv[3], "w") #output gene and pseu GFF file



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

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through pseu and make a dict of key: chromo_start_stop val: line and a dict
#of key: chromo val: [(start,stop),(start,stop),(start,stop)]
lineDict = {}
coordDict = {}
for line in pseu:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        chromo = lineLst[0]
        start, stop = OrderCoords(lineLst[3], lineLst[4])
        key = "%s_%s_%s" % (chromo,start,stop)
        lineDict[key] = line
        SaveIntoDict(chromo, (start,stop), coordDict)

#Add genes to both dicts where instead of lineDict having key: chrom_start_stop
#val: line for genes the val is val: nextLines (the line plus the mRNA, CDS and exon lines following)
lastLines = ""
key = ""
for line in gene:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "gene":
            chromo = lineLst[0]
            start, stop = OrderCoords(lineLst[3], lineLst[4])
            
            if not lastLines == "": #do nothing after first iteration
                lineDict[key] = lastLines
            
            key = "%s_%s_%s" % (chromo,start,stop)
            lastLines = line
            SaveIntoDict(chromo, (start,stop), coordDict)

        else:
            lastLines += line
#clean up and do the last line
else:
    lineDict[key] = lastLines

#make a chromoNumLst and sort it
chromoNumLst = []
for chromo in coordDict:
    chromoNum = chromo.strip("RrC")
    chromoNumLst.append(chromoNum)
chromoNumLst.sort()

#Go through chromoNumLst and coordDict, sort it and output everythin using lineDict
for chromoNum in chromoNumLst:
    chromo = "RrC" + str(chromoNum)
    coords = coordDict[chromo]
    coords.sort()
    for coord in coords:
        key = "%s_%s_%s" % (chromo, coord[0], coord[1])
        newLine = lineDict[key]
        out.write(newLine)





gene.close()
pseu.close()
out.close()
