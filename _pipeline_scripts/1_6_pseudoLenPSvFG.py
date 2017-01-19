#This script is designed to make a pseudogene comparative length file in the
#format:
#pseuName  pseuLen  funcGeneName  funcGeneLen  pseuLen/funcGeneLen
import sys

psSize = open(sys.argv[1])   #input 2col .size file of pseudogenes
fgSize = open(sys.argv[2])   #input 2col .size file of functional genes
out = open(sys.argv[3], "w") #output comparative length info file


def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)


        
#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through psSize and make a dict of key: fgName val: [(psName, size),]
psDict = {}
for line in psSize:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        psName = lineLst[0]
        size = lineLst[1]
        fgName = "_".join(psName.split("_")[1:])
        SaveIntoDict(fgName, (psName, size), psDict)

#Go through fgSize and output info
for line in fgSize:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        fgName = lineLst[0]
        fgLen = lineLst[1]
        try:
            for ps in psDict[fgName]:
                psName = ps[0]
                psLen = ps[1]
                sizeRatio = float(psLen)/float(fgLen)
                newLine = "%s\t%s\t%s\t%s\t%s\n" % (psName, psLen, fgName, fgLen, sizeRatio)
                out.write(newLine)
        except KeyError:
            continue




psSize.close()
fgSize.close()
out.close()
