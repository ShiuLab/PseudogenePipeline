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

#Takes group of names and makes all possible pairs from it (A,B,C --> A,B --> A,C --> B,C)
def GetPairs(words):
    #remove dashes from list
    BrGenes = []
    for word in words:
        if not word == "-":
            BrGenes.append(word)

    #make pairs
    if len(BrGenes) == 2:
        pairs = [(BrGenes[0], BrGenes[1]),]
    elif len(BrGenes) == 3:
        pairs = [(BrGenes[0], BrGenes[1]), (BrGenes[0], BrGenes[2]), (BrGenes[1], BrGenes[2])]
    else:
        print "error!!"
        
    return pairs

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
outSet = set()
for line in sub:
    pairsLst = []
    if not line.startswith("#"):
        lineLst = line.split("\t")
        subGenes = [lineLst[4], lineLst[5], lineLst[6].strip()]
        isValid = IsValid(subGenes)

        #if there are 2 or more genes for Br, continue processing
        if isValid == True:
            subPairs = GetPairs(subGenes)

            #b) go through orth and
            #    i) add whole lines to an output set if one of the orthologous
            #       generated pairs from 2a are represented in the line
            orth.seek(0)
            for ln in orth:
                if not ln.startswith("#"):
                    lnLst = ln.split("\t")
                    orthGenes = lnLst[2].split(",")
                    orthPairs = GetPairs(orthGenes)
                    orthPairs = MakeRevPairs(orthPairs)

                    #see if any subPairs are overlapping with orthPairs
                    for subPair in subPairs:
                        if subPair in orthPairs:
                            outSet.add(ln)
                    

#2) Go through ouput set and print the new lines
for line in outSet:
    out.write(line)





orth.close()
sub.close()
out.close()
