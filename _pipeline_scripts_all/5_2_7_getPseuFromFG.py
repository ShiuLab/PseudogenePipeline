#This script was designed to take a file with functional genes in it, a file
#with pseudogenes in it and make an output 1col file of the pseudogense
#associated with the functional genes
#Created by David E. Hufnagel on May 13, 2013
import sys

func = open(sys.argv[1])     #input finctional genes file (should have the genes of concern)
pseu = open(sys.argv[2])     #input psuedogenes file (should have all of them)
funCol = int(sys.argv[3])    #the column with the gene names (a col, not a python index)
psCol = int(sys.argv[4])     #the column with the pseudogene names (a col, not a python index)
out = open(sys.argv[5], "w") #output 1col file of pseudogenes


def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through pseu and make a dict of key: geneName val: pseuName
pseuDict = {}
for line in pseu:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        psName = lineLst[psCol-1]
	if psName.startswith("Ps"): #so it can handle a file with pseus and prots
            fgName = "_".join(psName.split("_")[1:])
	    SaveIntoDict(fgName, psName, pseuDict)

#Go through func and output pseudogene paralogs in the 1col format
for line in func:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        fgName = lineLst[funCol-1]
	if fgName in pseuDict:
	    for psName in pseuDict[fgName]:
                out.write(psName)
		out.write("\n")




func.close()
pseu.close()
out.close()
