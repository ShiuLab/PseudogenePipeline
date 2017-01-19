#this script is designed to take a file of stress related genes in the format:
#name up/down/no and a 4col of pseudogenes and make a 2col of stress related
#pseudogenes in the format:  name  up/down/no
#WARNING: SITUATION SPECIFIC

import sys

gene = open(sys.argv[1])     #input stress related genes file
pseu = open(sys.argv[2])     #input 4col pseudogenes file
out = open(sys.argv[3], "w") #output stress related pseudogene file


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through gene and make a dict of key: name val: up/down (no is excluded)
geneDict = {}
for line in gene:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        name = lineLst[0].replace(".","|")
        stress = lineLst[1]
        if not stress == "no":
            geneDict[name] = stress

#Go through pseu, find pseudogenes related to stress genes and output info
for line in pseu:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        psName = lineLst[0]
        gnName = psName.split("_")[1]
        if gnName in geneDict:
            newLine = "%s\t%s\n" % (psName, geneDict[gnName])
        else:
            newLine = "%s\tno\n" % (psName)
        out.write(newLine)




gene.close()
pseu.close()
out.close()
