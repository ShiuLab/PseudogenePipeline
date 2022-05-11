#This script is designed to filter out duplicate genes generated from tandem
#duplication from my 1:1:2:3 orthos file.  This is done by referencing the
#file from the Brassica paper with information about subgenomes and keeping
#only orthos with genes on multiple subgenomes.
#Created by David E. Hufnagel on July 13, 2012
"""Algorithm:
1) Go through sub and
    a) expand orthologous group into pairs(A,B,C --> A,B --> A,C --> B,C)
    b) go through orth and
        i) add whole lines to an output set if one of the orthologous
           generated pairs from 2a are represented in the line
2) Go through ouput set and print the new lines"""


import sys

orth = open(sys.argv[1])     #input orthos file
sub = open(sys.argv[2])      #subgenome info Br file
out = open(sys.argv[3], "w") #filtered output orthos file





#Returns true if the list is 2 or more genes long
def IsValid(BrGenes):
    if "-" in BrGenes:
        dashCnt = 0
        for gene in BrGenes:
            if gene == "-":
                dashCnt += 1
        if dashCnt <= 1:
            return True
        else:
            return False
    else:
        return True

#Despite the name it removes dashes
def GetSubGenes(words):
    geneLst = []
    for word in words:
        if word != "-":
            geneLst.append(word)

    return geneLst

#Takes a list of tuples (pairs) and returns a list with the reverse versions instead.
def MakeRevPairs(pairs):
    newPairs = []
    for pair in pairs:
        newPairs.append(pair)
        revPair = ((pair[1], pair[0]))
        newPairs.append(revPair)

    return newPairs





#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through sub and
#    a) expand orthologous group into pairs(A,B,C --> A,B --> A,C --> B,C)
badLst = []
for line in sub:
    pairsLst = []
    if not line.startswith("#"):
        lineLst = line.split("\t")
        subGenes = [lineLst[4], lineLst[5], lineLst[6].strip()]
        isValid = IsValid(subGenes)

        #if there are 2 or more genes in sub add them to the badLst
        if isValid == True:
            for gene in GetSubGenes(subGenes):
                badLst.append(gene)

for line in orth:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if not lineLst[2] in badLst:
            out.write(line)





orth.close()
sub.close()
out.close()
