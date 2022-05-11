#This script is designed to create a file for the Fisher Exact test in the
#format name\t1\t2\t3\t4 where 1=spe1+ 2=spe1- 3=spe2+ 4=spe2-
#Created by David E. Hufnagel on June 27, 2012

import sys

fish1 = open(sys.argv[1])    #input gene file for spe1 with GO terms for analysis
fish2 = open(sys.argv[2])    #input gene file for spe2 with GO terms for analysis
go = open(sys.argv[3])       #input GO tabular terms file
out = open(sys.argv[4], "w") #output fisher-exact file





def SaveIntoNumDict(key, dictX):
    if key not in dictX:
        dictX[key] = 1
    else:
        dictX[key] += 1

def ExtractFromFish(fd, special=""):
    speDict = {}
    geneSet = set()
    for line in fd:
        if not line.startswith("#"):
            if special == "At":
                #to exclude duplicated, chloroplast and mitochondrial genes
                if not line.startswith("AT") or line.startswith("ATC") or line.startswith("ATM"):
                    continue
            lineLst = line.split("\t")
            ID = lineLst[1].strip()
            if not ID.startswith("EC:"):  #causes a key error to do on the same line as the next line.
                if goDict[ID][1] == "biological_process" and goDict[ID][0] != "biological_process":
                    SaveIntoNumDict(ID, speDict)
##                elif goDict[ID][1] == "biological_process" and goDict[ID][0] == "biological_process":
##                    print "hello"
            geneSet.add(lineLst[0])
    geneNum = len(geneSet)
    return speDict, geneNum

def ReduceDict(dictx):
    newDict = {}
    for key in dictx:
        val = dictx[key]
        newVal = val[0]
        newDict[key] = newVal

    return newDict
        


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through go and put the file into a dict
goDict = {}    #a dict of key: GOID val: (name, namespace)
for line in go:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        goDict[lineLst[0]] = (lineLst[1], lineLst[2])#Go through go and put the file into a dict

#Go through fish1 and extract GO info for spe1 into a dict of key: name
#(obtained via the goDict) and number of times showing up in the file
#and the total number of genes.  Also reduce the key into just the name.
spe1Dict, spe1Tot = ExtractFromFish(fish1)

#The same deal for fish2
spe2Dict, spe2Tot = ExtractFromFish(fish2)#, "At")

#Remove namespaces from goDict
goDictR = ReduceDict(goDict)

#Go through the 2 dicts and output the information
for key1 in spe1Dict:
    val1 = spe1Dict[key1]
    if key1 in spe2Dict:
        val2 = spe2Dict[key1]
    else:
        val2 = "0"
    newLine = "%s\t%s\t%d\t%s\t%d\n" % \
              (goDictR[key1], val1, (int(spe1Tot)-int(val1)), val2, (int(spe2Tot)-int(val2)))
    out.write(newLine)



fish1.close()
fish2.close()
go.close()
out.close()
