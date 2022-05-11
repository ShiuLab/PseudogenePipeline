#This script was designed to take 2 orthologs files in the form At-Al-Br-Rr
#and make two 2col GO files: one for pseudogenes (dominant) and one for
#retained genes (2+ in Br and Rr AND no pseudogenes)
#Created by David E. Hufnagel on Sep 14, 2012
#WARNING FUNCTION CheckForPseudos IS SITUATION SPECIFIC FOR BR AND RR

import sys

pseu = open(sys.argv[1])        #input orthologs file that includes pseudogenes (singletons)
ret = open(sys.argv[2])         #input retained orthologs file
atGO = open(sys.argv[3])        #input At 2col GO file
rrPseuRef = open(sys.argv[4])   #input Rr pseudogenes 4 col file used as a reference for what functional genes go whith what pseudogenes
brPseuRef = open(sys.argv[5])   #input Br pseudogenes 4 col file used as a reference for what functional genes go whith what pseudogenes
psOut = open(sys.argv[6], "w")  #output 2col GO file for pseudogenized orthologous groups
retOut = open(sys.argv[7], "w") #output 2col GO file for retained orthologous groups

print
print "pseu orth:   %s" % (pseu)
print "ret orth:    %s" % (ret)
print "At GO:       %s" % (atGO)
print "Rr pseu 4col %s" % (rrPseuRef)
print "Br pseu 4col %s" % (brPseuRef)
print "pseu GO out: %s" % (psOut)
print "ret GO out:  %s\n" % (retOut)




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

#for when the pseudogenes need to be looked up in a dict of key: FG val: pseudo
def CheckForPseudos2(genes, dictx):
    hasPseudo = False
    for gene in genes:
        if gene in dictx:
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
retOut.write("#python %s\n" % (" ".join(sys.argv)))

#Go through atGO and import file into dict of key: AtName val: AtGO
goDict = {}
for line in atGO:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[0], lineLst[1].strip("\n"), goDict)

#Go through rrPseuRef and brPseuRef and import info into dict of key: FG name val: pseu name
rrBrPseuDict = {}
ImportPseuIntoDict(rrPseuRef, rrBrPseuDict)
ImportPseuIntoDict(brPseuRef, rrBrPseuDict)

#Go through pseu, filter out lines without any pseus in Br or Rr and put info
#into a list of At names (only 1 name per orthologous group)
pseuLst = []
for line in pseu:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        #Extract the names from the file
        brRrGenes = ExtractNames(lineLst)
        #Check if there are any pseudogenes in the file
        hasPseudo = CheckForPseudos2(brRrGenes, rrBrPseuDict)
        #if there are pseudogenes make the list of At gene names
        if not lineLst[0] == "[]" and hasPseudo == True:
            tempLst = ExtractAt(lineLst)
            pseuLst += tempLst
            #print line

#Go through ret, filter out lines with any pseus in Br or Rr and put into
#into a list of At names (only 1 name per orthologous group)
retLst = []
for line in ret:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        #Extract the names from the file
        brRrGenes = ExtractNames(lineLst)
        #Check if there are any pseudogenes in the file
        hasPseudo = CheckForPseudos2(brRrGenes, rrBrPseuDict)
        #if there aren't pseudogenes make the list of At gene names
        if not lineLst[0] == "[]" and hasPseudo == False:
            tempLst = ExtractAt(lineLst)
            retLst += tempLst

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
for gene in retLst:
    if gene in goDict:
        for go in goDict[gene]:
            newLine = "%s\t%s\n" % (gene, go)
            retOut.write(newLine)




pseu.close()
ret.close()
atGO.close()
rrPseuRef.close()
brPseuRef.close()
psOut.close()
retOut.close()
