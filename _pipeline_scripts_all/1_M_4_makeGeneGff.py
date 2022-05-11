#This script is designed to make a genic (not final) gff of Guarav's radish GFF
#containing only final genes
#Created by David E. Hufnagel on July 16, 2013
#NOTE: THIS IS THE VERSION OF STEP 1 THAT WORKS
#WARNING: HIGHLY SITUATION SPECIFIC
import sys

gff = open(sys.argv[1])      #input screwed up GFF file
four = open(sys.argv[2])     #input 4col of final genes
ref = open(sys.argv[3])      #input reference file for changing names
out = open(sys.argv[4], "w") #output file with only final genes (no alternative names)



#Go through four and make a list of good genes
goodLst = []
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        goodLst.append(lineLst[0])

#Go through ref and make a dict of key: altName val: realName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        altName = lineLst[0].split("-mRNA")[0]
        realName = lineLst[1]
        if realName in goodLst:
            refDict[altName] = lineLst[1]

#print len(refDict)

#Go through gff and output the proper lines
for line in gff:
    realName = "thisIsNotAName" #this helps realName reset in every iteration so the filter for good genes still works
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[1] == "maker":
            altName = lineLst[8].split("ID=")[1].split(";")[0].split("-mRNA")[0]
            if altName in refDict:
                realName = refDict[altName]
            if realName in goodLst:
                newLine = line.replace(altName,realName)
                out.write(newLine)




gff.close()
four.close()
ref.close()
out.close()
