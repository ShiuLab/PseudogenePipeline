#This script is desiged to take a parsed nodes file and determine from it how
#often Rr clusters with Br vs. clustering with Bn vs. being all alone.
#Created by David E. Hufnagel on oct 3, 2012
#WARNING: SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])      #input parsed nodes file
out = open(sys.argv[2], "w") #output file with relationship info




#does the work for left or right is neither is one gene long
def TestMultiGenes(genes, BrBn, BnRr, BrRr):
    speSet = set()
    for gene in genes:
        #make leftSpe
        if gene.startswith("rad") or gene.startswith("R") or gene.startswith("Rr"):
            speSet.add("Rr")
        elif gene.startswith("brs") or gene.startswith("B") or gene.startswith("Br"):
            speSet.add("Br")
        elif gene.startswith("bns") or gene.startswith("N") or gene.startswith("Bn"):
            speSet.add("Bn")
        else:
            print "\nerror here!\n"

    #test speSet to see if it's all the same spe
    if len(speSet) == 1:
        spe = list(speSet)[0]
        if spe == "Rr":
            BrBn += 1
        elif spe == "Br":
            BnRr += 1
        elif spe == "Bn":
            BrRr += 1
        else:
            print "\nerror here!\n"

    return BrBn, BnRr, BrRr





#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#set the counts for when two species are closer than another for a node
BrBn = 0
BrRr = 0
BnRr = 0
        
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")

        #get the left and right sides of the node
        left = lineLst[2].split(",")
        right = lineLst[3].split(",")

        #handle when one list is only 1 gene
        if len(left) == 1 or len(right) == 1:
            #set one
            if len(left) == 1:
                one = left[0]
            elif len(right) == 1:
                one = right[0]
                
            #use one to add to counts
            if one.startswith("rad") or one.startswith("R") or one.startswith("Rr"):
                BrBn += 1
            elif one.startswith("brs") or one.startswith("B") or one.startswith("Br"):
                BnRr += 1
            elif one.startswith("bns") or one.startswith("N") or one.startswith("Bn"):
                BrRr += 1
            else:
                print "\nerror here!\n"

        else:
            BrBn, BnRr, BrRr = TestMultiGenes(left, BrBn, BnRr, BrRr)
            BrBn, BnRr, BrRr = TestMultiGenes(right, BrBn, BnRr, BrRr)




out.write("BrBn: %s\n" % (BrBn))
out.write("BrRr: %s\n" % (BrRr))
out.write("BnRr: %s\n" % (BnRr))

        



inp.close()
out.close()
