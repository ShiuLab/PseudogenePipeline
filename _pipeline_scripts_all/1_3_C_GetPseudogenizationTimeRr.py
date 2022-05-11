#This script was designed to do the main and first calculations for pseudogene
#timing determination in Br-At orthologs.  More specifically it determines the
#timing of pseudogenization using the formula: 2*Fn*K(t-t1)+K*t1 = Nns where
#2*Fn*K(t-t1) is the nonsynonymous nucleotide substitutions built up before
#pseudogenization (1),  K*t1 is the nonysnonymous nucleotide substitutions
#built up after pseudogenization in the pseudogene (2) and Nns/site is the total
#nonsynonymous substitutions between the pseudogene and the closest related FG
#(3). Also: Fn*K = Ka between in and out. Fn = Ka/Ks between in and out. K =
#neutral evolution mutation rate. t = time since speciation. t1 = time since
#pseudogenization. (t-t1) = time before pseudogenization. Nns = num of
#nonsynonymous mutations per site.
#
#When organized for calculating t1, the formula goes:
#(2*Fn*K-Nns / K(Fn-1)) = t1
#
#This script creates a tab delimited info file in the format:
#PsName  InName  OutName  2*Fn*K(t-t1)  K*t1  Nns  t1  t
#Created by David E. Hufnagel on Dec 11, 2012.
#WARNING: FAIRLY SITUATION SPECIFIC

"""Tree for reference
         i
        1 i
       i   i
      2 i   i
     i-3-i   i
    Ps   in  out"""

import sys

inOut = open(sys.argv[1])    #input Ka/Ks info file for in vs. out gene
psIn = open(sys.argv[2])     #input Ka/Ks info file for pseudogene vs. in gene
psWGD = open(sys.argv[3])    #input info file having WGD derived pseudogenes in the first col
orthos = open(sys.argv[4])   #input orthologs file to get Br-At relationships 
out = open(sys.argv[5], "w") #output info file
K = float(sys.argv[6])       #pseudogene evolution rate (usually assumed to be neutral .007)
###WARNING: CHECK STEP2###



#In an orthos file, extract the genes for one col (one spe)
#The assumed format is: [At1, At2] [Al1,Al2] [Br1, Br2] [Rr1,Rr2]]
def ExtractGenesFromCol(chunk):
        genes = chunk[2:-2].split(",")
        genes2 = []
        for gene in genes:
            genes2.append(gene.strip(" ").strip("'"))

        return genes2

#Given two lists, creates all the possible combinatorial pairs.  A = [1,2,3]
#B = [6,7] GetAllPairs(A,B) --> [[1,6],[1,7],[2,6],[2,7],[3,6],][3,7]
def GetAllCombs(listA,listB):
    newLst = []
    for geneA in listA:
        for geneB in listB:
            newLst.append([geneA, geneB])

    return newLst

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)
            



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
out.write("#  PsName  InName  OutName  2*Fn*K(t-t1)  K*t1  Nns  t1  t t1/t\n")

#1) Put Rr-At orthologs from orthos into orthoDict. key: RrGene val: (AtGene1, AtGene2)
#assumes orthos file is in the form: [At1, At2] [Al1,Al2] [Br1, Br2] [Rr1,Rr2]]
orthoDict = {}
cnt = 0
for line in orthos:
    if not line.startswith("#"):
        #extract At and Rr genes
        lineLst = line.strip().split("\t")
        atGenes = ExtractGenesFromCol(lineLst[0])
        rrGenes = ExtractGenesFromCol(lineLst[3])

        if rrGenes != [''] and atGenes != ['']:
            #make all possible gene pairs between these two groups
            pairs = GetAllCombs(rrGenes, atGenes)

            #put info into orthoDict
            for pair in pairs:
                SaveIntoDict(pair[0], pair[1], orthoDict)

#2) Go through psIn and...
psDict = {} #dict of key: psName val: [[RrGene, AtGene, PsRr:ka, PsRr:ks, PsRr:ka/ks],]
rrDict = {} #dict of key: RrGene val: psName
for line in psIn:
    if not line.startswith("#"):
        #get pairs
        lineLst = line.strip().split("\t")
        ### MAKE SURE THESE TWO NAMES ARE SET RIGHT ###
        rrName = lineLst[1]
        psName = lineLst[0]
        
        if rrName.startswith("RrC") and psName.startswith("PsRr") and rrName in orthoDict:
            atNames = orthoDict[rrName]
            for atName in atNames:        
                #a) Get pseudos related to Rr gene and the pair's ka, ks and
                #   ka/ks.
                ####This is what makes _mod different ####
                psRrKa = lineLst[4]; psRrKs = lineLst[5]; psRrKaKs = lineLst[6]

                #b) Put info into PsDict with key: psName val: [[RrGene, AtGene, PsRr:ka,
                #   PsRr:ks, PsRr:ka/ks],]
                SaveIntoDict(psName, [rrName, atName, psRrKa, psRrKs, psRrKaKs], psDict)

            #c) Make BrDict of key: RrGene val: psName
            SaveIntoDict(rrName, psName, rrDict)

#3) Go through inOut and make inOutDict with key: rrName;AtName val: (Ka, Ks, Ka/Ks) and make divDict with key: rrName;AtName val: divTime
inOutDict = {}
divDict = {}
for line in inOut:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        name1 = lineLst[0]
        name2 = lineLst[1]
        if (name1.startswith("AT") or name2.startswith("AT")) and (name1.startswith("RrC") or name2.startswith("RrC")):
            if name1.startswith("AT"):
                atName = name1
                rrName = name2
            elif name2.startswith("AT"):
                atName = name2
                rrName = name1

            key = "%s;%s" % (rrName, atName)
            inOutDict[key] = (lineLst[4], lineLst[5], lineLst[6])

            #Make divDict using Ks values to determine species divergence: T = Ks/2*K 
            divTime = float(lineLst[5]) / (2.0*K)
            divDict[key] = (divTime)

            

#4) Go through psWGD and make psList to contain WGD derived pseudogenes
psList = []
for line in psWGD:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        rrName = "_".join(lineLst[0].split("_")[1:])
        psList.append(rrName)

#4) Go through psDict
#a) Calculate Fn*K(t-t1), K*t1 and t1
#b) output info in the format:
#   PsName  InName  OutName  Fn*K(t-t1)  K*t1  Nns  t1  t
for psName in psDict.keys():
    for group in psDict[psName]:
        rrName = group[0]
        atName = group[1]

        #a) Calculate Fn*K(t-t1), K*t1 and t1
        key = "%s;%s" % (rrName, atName)
        if key in inOutDict and key in divDict and rrName in psList:
            Fn = float(inOutDict[key][2]) #Ka/Ks between in and out
            Nns = float(group[2])         #Nns = num of nonsynonymous mutations per site between in and Ps 
            t = float(divDict[key])
            t1 = ((2*Fn*K*t)-Nns) / K*(Fn-1.0)
            ratio = t1/t
            one = 2*Fn*K*(t-t1)
            two = K*t1
            
            #b) output info in the format:
            #  PsName  InName  OutName  2*Fn*K(t-t1)  K*t1  Nns  t1  t t1/t
            newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                      (psName, rrName, atName, one, two, Nns, t1, t, ratio)
            out.write(newLine)
        



inOut.close()
psIn.close()
psWGD.close()
orthos.close()
out.close()

