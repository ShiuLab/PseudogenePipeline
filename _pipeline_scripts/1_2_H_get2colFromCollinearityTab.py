#This script is designed to make a 2col homology-like file from the tabular
#collinearity file with Gaurav's style
#Created by David E. Hufnagel on Jan 23, 2013

import sys

coll = open(sys.argv[1])     #input tabular Gaurav-style collinearity file 
four = open(sys.argv[2])     #input 4col genes file (must contain all genes in blocks)
out = open(sys.argv[3], "w") #output 2col file



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through 4col file and make refDict of key = geneName val = chromo
refDict = {}
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        refDict[lineLst[0]] = lineLst[1]

#Go through collinearity and make a set of newLines
newLineSet = set()
for line in coll:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        blockID = lineLst[2].split("|")[0]
        chromo1 = refDict[lineLst[0]]
        chromo2 = refDict[lineLst[1]]
        newLine = "%s_block_%s_%s|%s\t%s_block_%s_%s|%s\n" % (chromo1, blockID, chromo1, chromo2, chromo2, blockID, chromo1, chromo2)
        newLineSet.add(newLine)

#Output info
for newLine in newLineSet:
    out.write(newLine)
    



coll.close()
four.close()
out.close()
