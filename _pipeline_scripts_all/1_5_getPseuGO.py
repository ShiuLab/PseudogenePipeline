#This script is designed to make a 4 column file of pseudogenes and the GO
#categories of their functional counterparts in GOID form
#Created by David E. Hufnagel on June 28, 2012
#WARNING: WILL NEED ADJUSTMENT IF A NEW SPECIES IS USED
"""Algorithm:
1) Go through go and save the file into a goDict
2) Go through goRef and save the file info into a dict.
3) Go through pseu and use the goDict to put the data together and output
the info into out"""

import sys

pseu = open(sys.argv[1])     #input 4col pseudogene file with final code-derived names    (.real4col.RMfilt.cdnm)
go = open(sys.argv[2])       #input 2col GO file for functional genes                     (?)
goRef = open(sys.argv[3])    #input 2col GO reference file with IDs, names and namespaces (?)
out = open(sys.argv[4], "w") #output pseudogene GO file                                   (.go)
spe = sys.argv[5]            #the 2 character species ID (for At)                         (Ex: At)





def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)




        
#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
#1) Go through go and save the file into a dict. key: protName val: GOID
goDict = {}
for line in go:
    if not line.startswith("#"):
        if spe == "At":
            if not line.startswith("AT") or line.startswith("ATC") or line.startswith("ATM"):
                continue
        lineLst = line.split("\t")
        if not lineLst[1].startswith("EC:"):  #for EC annotations (Rr)
            prot = lineLst[0].split("|")[0]   #also for Rr
            SaveIntoDict(prot, lineLst[1].strip(), goDict)

#print len(goDict)

#2) Go through goRef and save the file info into a dict. key: GOID val: (name, namespace)
refDict = {}
for line in goRef:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        #this says if the name is a namespace (lazy annotaters!) don't bother adding it to the dict
        if not lineLst[1] == "molecular_function" and\
           not lineLst[1] == "biological_process" and\
           not lineLst[1] == "cellular_component":
            refDict[lineLst[0]] = (lineLst[1], lineLst[2])

#print len(refDict)

#3) Go through pseu and use the goDict and refDict to put the data together
#and output the info into out
for line in pseu:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        prot = "_".join(lineLst[0].split("_")[1:])
        if prot in goDict:
            pseud = lineLst[0]
            for GOID in goDict[prot]:
                #sometimes a GOID is not in the refDict because I didn't import things there the name is a namespace
                try:
                    GOname = refDict[GOID][0]
                except:
                    continue
                namespace = refDict[GOID][1]
                newLine = "%s\t%s\t%s\t%s\n" % (pseud, GOID, GOname, namespace)
                out.write(newLine)




pseu.close()
go.close()
out.close()
