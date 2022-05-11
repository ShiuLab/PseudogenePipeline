#This script is designed to filter out duplicate genes generated from tandem
#duplication from my orthos files.  This is done by referencing the
#file from the Brassica paper with information about subgenomes and keeping
#only orthos with genes in the file and making 2 kinds of outputs 1) the orth
#lines where one or more genes is present in sub and 2) the orth lines where
#all genes were present and in the same group in sub.
#Created by David E. Hufnagel on July 16, 2012
"""Algorithm:
1) Go through sub and...
   a) Make a list of genes in sub
   b) Make a list of lists with the same genes in groups based on what line
      they were found on (groups of 1-3)
2) Go through orth and...
   a) Output lines with genes from the subLst into outMain
   b) Make a list of lists of orth genes.
3) Go through orthLst and output lines where genes in orthLst are in the same
   groups in subLst (this includes making pairs and reverse pairs and comparing
   them)"""


import sys

orth = open(sys.argv[1])         #input orthos file
sub = open(sys.argv[2])          #subgenome info Br file
outMain = open(sys.argv[3], "w") #filtered output orthos file
outSame = open(sys.argv[4], "w") #filtered output orthos file where orthos were the same in both input files ("the same" means all genes from orth)





#Takes group of names and makes all possible pairs from it (A,B,C --> A,B --> A,C --> B,C)
def GetPairs(words):
    #remove dashes from list
    BrGenes = []
    for word in words:
        if not word == "-":
            BrGenes.append(word)

    #make pairs
    if len(BrGenes) == 1:
        pairs = [(BrGenes[0],""),]
    elif len(BrGenes) == 2:
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

def RemoveDashes(genes):
    outLst = []
    for gene in genes:
        if gene != "-":
            outLst.append(gene)

    return outLst
        





#Write the users command line prompt on the first line of the output files.
outMain.write("#python %s\n" % (" ".join(sys.argv)))
outSame.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through sub and...
#   a) Make a list of genes in sub
#   b) Make a list of lists with the same genes in groups based on what line
allSubGenes = []
subGroups = []
for line in sub:
    if not line.startswith("#"):
        lineLst = line.split("\t")

        #get genes and group of genes if lines are not empty
        if not (lineLst[4] == "-" and lineLst[5] == "-" and lineLst[6].strip() == "-"):
            subGenes = [lineLst[4], lineLst[5], lineLst[6].strip()]
            subGenes = RemoveDashes(subGenes)
            subGroups.append(subGenes)

            for gene in subGenes:
                allSubGenes.append(gene)

#2) Go through orth and...
#   a) Output lines with genes from the subLst into outMain
#   b) Make a list of lists of orth genes.
#   c) Make a dict of key: orthGroup val: line
orthGroups = []
orthDict = {}
for line in orth:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        orthGenes = lineLst[2].split(",")

        #a) Output lines with genes from the subLst into outMain
        doWrite = False
        for orthGene in orthGenes:
            if orthGene in allSubGenes:
                doWrite = True
        if doWrite == True:
            outMain.write(line)

        #b) Make a list of lists of orth genes.
        orthGroups.append(orthGenes)

        #c) Make a dict of key: orthGroup val: line
        orthDict[str(orthGenes)] = line

#3) Go through subLst and make subPairs from sunGroups
#   including reverse pairs (A,B,C --> A,B  A,C  B,C  B,A  C,A  C,B)
allSubPairs = []
for subGroup in subGroups:
    subFPairs = GetPairs(subGroup)  #get forward pairs
    subPairs = MakeRevPairs(subFPairs) #add reverse pairs
    allSubPairs += subPairs

#4) Go through orthGroups and output lines where genes in orthLst are in the same
#   groups as in subLst (includes making pairs for orthGroups)
for orthGroup in orthGroups:
    orthFPairs = GetPairs(orthGroup)  #get forward pairs (reverse not needed becasue it's already done in sub)

    doWrite = False
    for orthPair in orthFPairs:
        if orthPair in allSubPairs:
            doWrite = True
    if doWrite == True:
        newLine = orthDict[str(orthGroup)]
        outSame.write(newLine)
            


    

        





orth.close()
sub.close()
outMain.close()
outSame.close()
