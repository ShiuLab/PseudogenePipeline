#This script was designed to take 2 orthologs files in the form At-Al-Br-Rr
#and make two 2col GO files: one for pseudogenes (dominant) and one for
#retained genes (2+ in Br and Rr AND no pseudogenes)
#Created by David E. Hufnagel on Sep 14, 2012
#WARNING FUNCTION CheckForPseudos IS SITUATION SPECIFIC FOR BR AND RR

import sys

pseu = open(sys.argv[1])         #input orthologs file that includes pseudogenes (singletons)
atGO = open(sys.argv[2])         #input At 2col GO file
psOut = open(sys.argv[3], "w")   #output 2col GO file for pseudogenized orthologous groups
noPsOut = open(sys.argv[4], "w") #output 2col GO file for retained orthologous groups

print
print "pseu orth:      %s" % (pseu)
print "At GO:          %s" % (atGO)
print "pseu GO out:    %s" % (psOut)
print "noPseu GO out:  %s\n" % (noPsOut)




def ExtractAt(lineLst):
    atLst = []
    at = lineLst[0].strip("[]").split(",")
    for gene in at:
        gene = gene.strip(" ''")
        atLst.append(gene)

    return atLst

def ExtractNames(lineLst):
    brRrLst = []
    br = lineLst[2].strip("[]").split(",")
    for gene in br:
        gene = gene.strip(" ''")
        brRrLst.append(gene)

    rr = lineLst[3].strip("[]\n").split(",")
    for gene in rr:
        gene = gene.strip(" ''")
        brRrLst.append(gene)
        
    return brRrLst

#for when pseudogenes would be included in the list
def CheckForPseudos(genes):
    hasPseudo = False
    for gene in genes:
        if "PsBr" in gene or "PsRr" in gene:
            hasPseudo = True

    return hasPseudo

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

def ImportPseuIntoDict(fd, dictx):
    for line in fd:
        if not line.startswith("#"):
            lineLst = line.split("\t")
            pseuName = lineLst[0].split("_")[0]
            fGName = "_".join(lineLst[0].split("_")[1:])
            SaveIntoDict(fGName, pseuName, dictx)
            

        

#Write the users command line prompt on the first line of the output files.
psOut.write("#python %s\n" % (" ".join(sys.argv)))
noPsOut.write("#python %s\n" % (" ".join(sys.argv)))

#Go through atGO and import file into dict of key: AtName val: AtGO
goDict = {}
for line in atGO:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[0], lineLst[1].strip("\n"), goDict)

#Go through pseu, filter out lines without any pseus in Br or Rr and put info
#into a list of At names (only 1 name per orthologous group)
pseuLst = []
noPsLst = []
for line in pseu:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        #Extract the names from the file
        brRrGenes = ExtractNames(lineLst)
        #Check if there are any pseudogenes in the file
        hasPseudo = CheckForPseudos(brRrGenes)
        #if there are pseudogenes make the list of At gene names
        if not lineLst[0] == "[]" and hasPseudo == True:
            tempLst = ExtractAt(lineLst)
            pseuLst += tempLst
        elif not lineLst[0] == "[]" and hasPseudo == False:
            tempLst = ExtractAt(lineLst)
            noPsLst += tempLst

#Go through pseuLst and output gene names with GO categories into 2col GO output file
for gene in pseuLst:
    if gene.startswith("Ps"):
        key ="_".join(gene.split("_")[1:])
    else:
        key = gene
    if key in goDict:
        for go in goDict[key]:
            newLine = "%s\t%s\n" % (gene, go)
            psOut.write(newLine)

#Go through retLst and output gene names with GO categories into 2col GO output file
for gene in noPsLst:
    if gene in goDict:
        for go in goDict[gene]:
            newLine = "%s\t%s\n" % (gene, go)
            noPsOut.write(newLine)




pseu.close()
atGO.close()
psOut.close()
noPsOut.close()
