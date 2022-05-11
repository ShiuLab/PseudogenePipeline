#This script is designed to take my differential expression files generated from
#parsed edgeR output files (format: ##info lines followed by a list of \n
#delimited gene names) and add info from a few different sources to make the
#format:
#RrGeneName   AtOrthoGeneName   GO_ID   GO type   GO description
#Created by David E. Hufnagel on Nov 25, 2012
"""Algorithm:
1) Go through rrGenes and make a list of all gene names
2) Go through orthos and make a dict of key: atGene val: rrGene for all rrGenes
   in rrList
3) Go through atGOs and through atDict and add GOID to the dict with GOID
   as the new key
4) Go through goTab and through atDict and add GO type and GO description to
   the dict info in the form of a list of newLines
5) Go though atDict and rrGenes (again) and output info in the format:
RrGeneName   AtOrthoGeneName   GO_ID   GO type   GO description"""
import sys

rrGenes = open(sys.argv[1])  #EdgeR parsed output file for Rr
orthos = open(sys.argv[2])   #pairwise 2col orthologs file
atGOs = open(sys.argv[3])    #whole genome At 2col GO file
goTab = open(sys.argv[4])    #tabular GO ID database file
out = open(sys.argv[5], "w") #output 5col info file




#1) Go through rrGenes and make a list of all gene names
rrSet = set([])
for line in rrGenes:
    if not line.startswith("#"):
        rrSet.add(line.strip())
rrList = list(rrSet)
        
#2) Go through orthos and make a dict of key: atGene val: rrGene for all rrGenes
#   in rrList
atDict = {}
for line in orthos:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        #work with only legitimite lines
        if lineLst[0] in rrList and lineLst[1].startswith("AT"):
            atDict[lineLst[1]] = lineLst[0]
        elif lineLst[1] in rrList and lineLst[0].startswith("AT"):
            atDict[lineLst[0]] = lineLst[1]

#3) Go through atGOs and through atDict and add GOID to the dict with GOID
#   as the new key
goIDDict = {}
for line in atGOs:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[0] in atDict:
            if lineLst[1] not in goIDDict:
                goIDDict[lineLst[1]] = [(lineLst[0], atDict[lineLst[0]])]
            else:
                goIDDict[lineLst[1]].append((lineLst[0], atDict[lineLst[0]]))

#4) if there is no GOID for the particular At ortho, add that to the goIDDict
for atGene in atDict:
    doAdd = True
    for go in goIDDict.values():
        for group in go:
            if atGene == group[0]:
                doAdd = False
            
    if doAdd == True:
        #print atGene
        if "-" not in goIDDict:
            goIDDict["-"] = [(atGene, atDict[atGene])]
        else:
            goIDDict["-"].append((atGene, atDict[atGene]))

#5) Go through goTab and through atDict and add GO type and GO description to
#   the dict info in the form of a list of newLines
outLst = []
for line in goTab:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[0] in goIDDict:
            for group in goIDDict[lineLst[0]]:
                rrGene = group[1]
                atGene = group[0]
                goID = lineLst[0]
                goType = lineLst[2]
                goDesc = lineLst[1] #GO description
                newLine = "%s\t%s\t%s\t%s\t%s\n" % (rrGene, atGene, goID, goType, goDesc)
                outLst.append(newLine)

#6) Add lines without GOIDs to the outLst
for group in goIDDict["-"]:
    rrGene = group[1]
    atGene = group[0]
    goID = "-"
    goType = "-"
    goDesc = "-" #GO description
    newLine = "%s\t%s\t%s\t%s\t%s\n" % (rrGene, atGene, goID, goType, goDesc)
    outLst.append(newLine)
                

#7) Go though atDict and rrGenes (again) and output info in the format:
#RrGeneName   AtOrthoGeneName   GO_ID   GO type   GO description
rrGenes.seek(0)
for line in rrGenes:
    if not line.startswith("#"):
        isDone = False
        for newLine in outLst:
            if line.strip() in newLine:
                out.write(newLine)
                isDone = True
        #Handle when there is no At counterpart to the Rr gene
        if isDone == False:
            out.write("%s\t-\t-\t-\t-\n" % (line.strip()))
    else:
        out.write(line)
            




rrGenes.close()
orthos.close()
atGOs.close()
goTab.close()
out.close()
