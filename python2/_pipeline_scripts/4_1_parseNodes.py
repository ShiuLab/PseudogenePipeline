#This script is desiged to parse a nodes file (a file derived from a newick
#phylogeny file using 4_1_makeNodesTabFromNewick.py) to keep only nodes with
#3 or more genes, all three species represented (Br,Bn,Rr) and never the same
#species on both sides (left and right)
#Created by David E. Hufnagel on Oct 3, 2012
#WARNING: SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])      #the input unparsed nodes file
out = open(sys.argv[2], "w") #the output parsed nodes file




def GetSpe(genes):
    spe = set()
    for gene in genes:
        if gene.startswith("rad") or gene.startswith("R") or gene.startswith("Rr"):
            spe.add("Rr")
        elif gene.startswith("brs") or gene.startswith("B") or gene.startswith("Br"):
            spe.add("Br")
        elif gene.startswith("bns") or gene.startswith("N") or gene.startswith("Bn"):
            spe.add("Bn")
        else:
            print "\nerror here!\n"

    return list(spe)




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        
        #determine whether to keep and therefore to write the line 
        doWrite = True
        speLst = []
        for spe in lineLst[1].split(","):
            speLst.append(spe)

        if len(speLst) < 3:  ##this test that all species are there AND that there are 3+ genes simultaneously
            doWrite = False

        ##test whether any species is represented in left and right  
        leftGenes = lineLst[2].split(",")
        rightGenes = lineLst[3].split(",")
        rightGenes[-1] = rightGenes[-1].strip()

        leftSpe = GetSpe(leftGenes)
        rightSpe = GetSpe(rightGenes)

        testLst = []
        for spe in leftSpe:
            if spe in rightSpe:
                doWrite = False

        #output info
        if doWrite == True:
            out.write(line)




inp.close()
out.close()
