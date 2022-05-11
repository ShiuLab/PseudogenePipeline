#This script is designed to separate pseudogenes with stop codons or frameshifts
#from those without and to output the info in the form of two 4col files

import sys

inp = open(sys.argv[1])          #input 4col pseudogene file with full names (.real4col.RMfilt)
ref = open(sys.argv[2])          #input reference file with full names and code names
outStop = open(sys.argv[3], "w") #output 4col file of pseudogenes with stop codons or frameshifts
outNone = open(sys.argv[4], "w") #output 4col file of pseudogenes without stop codons and frameshifts




#writes the users command line prompt on the first line of the output file.
outStop.write("#python %s\n" % (" ".join(sys.argv)))
outNone.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and get dict of key: longPsName val: codeName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        refDict[lineLst[1]] = lineLst[0]

#Go through inp and determine which pseudogenes belong in which category and
#output the info into two lists with code names
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        codeLst = lineLst[0].split(";")[3].split("|")
        codeName = refDict[lineLst[0]]
        newLine = "%s\t%s\t%s\t%s\n" % (codeName, lineLst[1], lineLst[2], lineLst[3])
        #check if it is in the none group
        if codeLst[0] == "0" and codeLst[1] == "0" and codeLst[2] == "0" and codeLst[3] == "0":
           outNone.write(newLine)
        #if it isn't in the none group
        else:
           outStop.write(newLine)





inp.close()
outStop.close()
outNone.close()
