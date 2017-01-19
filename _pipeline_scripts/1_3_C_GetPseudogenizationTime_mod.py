#This script was designed to do the main and first calculations for pseudogene
#timing determination in Br-At orthologs.  More specifically it determines the
#timing of pseudogenization using the formula: Fn*K(t-t1)+K*t1 = Nns where
#Fn*K(t-t1) is the nonsynonymous nucleotide substitutions built up before
#pseudogenization (1),  K*t1 is the nonysnonymous nucleotide substitutions
#built up after pseudogenization in the pseudogene (2) and Nns/site is the total
#nonsynonymous substitutions between the pseudogene and the closest related FG
#(3). Also: Fn*K = Ka between in and out. Fn = Ka/Ks between in and out. K =
#neutral evolution mutation rate. t = time since speciation. t1 = time since
#pseudogenization. (t-t1) = time before pseudogenization. Nns = num of
#nonsynonymous mutations per site.
#
#When organized for calculating t1, the formula goes:
#(Fn*K-Nns / K(Fn-1)) = t1
#
#This script creates a tab delimited info file in the format:
#PsName  InName  OutName  Fn*K(t-t1)  K*t1  Nns  t1  t t1/t  fnkt-Nns  K*(Fn-1) (Fn,Nns,Fnkt)
#Created by David E. Hufnagel on Dec 11, 2012.
#The _mod version was created on Jan 20, 2013 to use different input Ka/Ks
#files.  The files can be found in the HPCC folder /mnt/home/hufnag30/Shiu/
#1_RadishProject/3_AtPseudogeneTiming/C_ActualTiming/1_ChouMethod/2_Jan17_2013
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
out.write("#  PsName  InName  OutName  Fn*K(t-t1)  K*t1  Nns  t1  t t1/t  fnkt-Nns  K*(Fn-1)  (Fn,Nns,Fnkt)\n")

#1) Put Br-At orthologs from orthos into orthoDict. key: BrGene val: (AtGene1, AtGene2)
#assumes orthos file is in the form: [At1, At2] [Al1,Al2] [Br1, Br2] [Rr1,Rr2]]
orthoDict = {}
cnt = 0
for line in orthos:
    if not line.startswith("#"):
        #extract At and Br genes
        lineLst = line.strip().split("\t")
        atGenes = ExtractGenesFromCol(lineLst[0])
        brGenes = ExtractGenesFromCol(lineLst[2])

        if brGenes != [''] and atGenes != ['']:
            #make all possible gene pairs between these two groups
            pairs = GetAllCombs(brGenes, atGenes)

            #put info into 
            for pair in pairs:
                SaveIntoDict(pair[0], pair[1], orthoDict)

#2) Go through psIn and...
psDict = {} #dict of key: psName val: [[BrGene, AtGene, PsBr:ka, PsBr:ks, PsBr:ka/ks],]
brDict = {} #dict of key: BrGene val: psName
for line in psIn:
    if not line.startswith("#"):
        #get pairs
        lineLst = line.strip().split("\t")
        brName = lineLst[0]
        psName = lineLst[1]
        
        if brName.startswith("Bra") and psName.startswith("PsBr") and brName in orthoDict:
            atNames = orthoDict[brName]
            for atName in atNames:        
                #a) Get pseudos related to Br gene and the pair's ka, ks and
                #   ka/ks.
                ####This is what makes _mod different ####
                psBrKa = lineLst[2]; psBrKs = lineLst[3]
                
                #make sure to handle division by 0
                if float(psBrKs) == 0:
                    psBrKaKs = (float(psBrKa) / 0.000000001)
                else:
                    psBrKaKs = (float(psBrKa) / float(psBrKs))

                #b) Put info into PsDict with key: psName val: [[BrGene, AtGene, PsBr:ka,
                #   PsBr:ks, PsBr:ka/ks],]
                SaveIntoDict(psName, [brName, atName, psBrKa, psBrKs, psBrKaKs], psDict)

            #c) Make BrDict of key: BrGene val: psName
            SaveIntoDict(brName, psName, brDict)

#3) Go through inOut and make inOutDict with key: brName;AtName val: (Ka, Ks, Ka/Ks) and make divDict with key: brName;AtName val: divTime
inOutDict = {}
divDict = {}
for line in inOut:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        name1 = lineLst[0]
        name2 = lineLst[1]
        if (name1.startswith("AT") or name2.startswith("AT")) and (name1.startswith("Bra") or name2.startswith("Bra")):
            if name1.startswith("AT"):
                atName = name1
                brName = name2
            elif name2.startswith("AT"):
                atName = name2
                brName = name1

            key = "%s;%s" % (brName, atName)
            inOutDict[key] = (lineLst[4], lineLst[5], lineLst[6])

            #Make divDict using Ks values to determine species divergence: T = Ks/2*K 
            divTime = float(lineLst[5]) / (2.0*K)
            divDict[key] = (divTime)

            

#4) Go through psWGD and make psList to contain WGD derived pseudogenes
psList = []
for line in psWGD:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        brName = lineLst[0].split("_")[1]
        psList.append(brName)

#4) Go through psDict
#a) Calculate Fn*K(t-t1), K*t1 and t1
#b) output info in the format:
#   PsName  InName  OutName  Fn*K(t-t1)  K*t1  Nns  t1  t
for psName in psDict.keys():
    for group in psDict[psName]:
        brName = group[0]
        atName = group[1]

        #a) Calculate Fn*K(t-t1), K*t1 and t1
        key = "%s;%s" % (brName, atName)
        if key in inOutDict and key in divDict and brName in psList:
            Fn = float(inOutDict[key][2]) #Ka/Ks between in and out
            Nns = float(group[2])         #Nns = num of nonsynonymous mutations per site between in and Ps 
            t = float(divDict[key])
            t1 = ((Fn*K)-Nns) / K*(Fn-1.0)
            ratio = t1/t
            one = Fn*K*(t-t1)
            two = K*t1

            fnkt = Fn*K*t
            num1 = fnkt-Nns
            den1 = K*(Fn-1)
            str12 = "(%s,%s,%s)" % (Fn,Nns,fnkt)
            
            #b) output info in the format:
            #  PsName  InName  OutName  Fn*K(t-t1)  K*t1  Nns  t1  t t1/t  fnkt-Nns  K*(Fn-1) (Fn,Nns,Fnkt)
            newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                      (psName, brName, atName, one, two, Nns, t1, t, ratio, num1, den1, str12)
            out.write(newLine)
        



inOut.close()
psIn.close()
psWGD.close()
orthos.close()
out.close()

