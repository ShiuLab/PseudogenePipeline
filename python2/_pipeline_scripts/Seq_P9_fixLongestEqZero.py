#This script is designed to fix the "longest=0 problem"  which is where a gff
#file has genes with only one mRNA, and the mRNA is labeled as longest=0

import sys

gff = open(sys.argv[1])
out = open(sys.argv[2], "w")



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through gff and make a dict of key: geneName val: list of mRNA info lines
gffDict = {}
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "gene":
            if gffDict != {}:
                lastGene = currGene
                gffDict[lastGene] = lastMRNAs
            else:
                gffDict["deleteMe"] = 0 #so that the if gffDict != 0 statement becomes false after the first iteration
            lastMRNAs = []
            currGene = lineLst[8].split("Name=")[1]
        if lineLst[2] == "mRNA":
            lastMRNAs.append(lineLst[8])
#get last line
else:
    lastGene = currGene
    gffDict[lastGene] = lastMRNAs


#remove false item
gffDict.pop("deleteMe")

#Go through dict and if len(val) is 1 and longest=0 is contained in the info line, save the mRNA name in a toChangeLst
toChangeLst = []
for geneName in gffDict:
    mRNAinfo = gffDict[geneName]
    if len(mRNAinfo) == 1:
        if "longest=0" in mRNAinfo[0]:
            mRNAname = mRNAinfo[0].split("Name=")[1].split(";")[0]
            toChangeLst.append(mRNAname)
        
    elif len(mRNAinfo) == 0:
        print "ERROR, NO MRNA'S FOR THIS GENE:", geneName
        
#    else:
#        print geneName

#print toChangeLst

#Go through gff again and change longest=0 to longest=1 on mRNA lists with names in the toChangeLst
gff.seek(0)
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "mRNA":
            currName = lineLst[8].split("Name=")[1].split(";")[0]
            if currName in toChangeLst:
                line = line.replace("longest=0","longest=1")

        out.write(line)


gff.close()
out.close()
