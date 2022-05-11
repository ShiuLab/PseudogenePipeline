#This script is designed to go through Guarav's radish gff file and get rid of
#code names in a thoughtful manor
#Created by David E. Hufnagel on July 3, 2013
#WARNING: HIGHLY SITUATION SPECIFIC
import sys

inp = open(sys.argv[1])      #input gff file with old names
ref = open(sys.argv[2])      #reference file with old and new names
out = open(sys.argv[3], "w") #output gff file with new names


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and make a dict of key: oldName val: newName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        refDict[lineLst[0]] = lineLst[1]

#Go through inp
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        
        #get oldName and determine nweName based on refDict and lineLst
        oldName = lineLst[8].split("ID=")[1].split(":")[0] #for exon, mRNA and CDS
        if lineLst[2] == "gene":
            oldName = lineLst[8].split("ID=")[1].split(";")[0] #not in ref
            for refName in refDict:
                if oldName in refName:
                    newName = refDict[refName]

        elif lineLst[2] == "exon":
            endPiece =  ":" + "-".join(lineLst[8].split("ID=")[1].split(";")[0].split("-")[-2:])
            newName = refDict[refName] + endPiece
        elif lineLst[2] == "mRNA":
            endPiece =  ":" + lineLst[8].split("ID=")[1].split(";")[0].split("-")[-1]
            newName = refDict[refName] + endPiece
            print lineLst
            print oldName
            print newName
            print
        elif lineLst[2] == "CDS":
            pass

        
        








inp.close()
ref.close()
out.close()
