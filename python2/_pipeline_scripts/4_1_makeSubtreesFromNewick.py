#This script is designed to make a subtree from a Newick phylogeny file.
#The output file is all nodes with 3 or more genes, all three species
#represented (Br,Bn,Rr) and never the same species on both sides (left and
#right).  This script is largely derived from 4_1_makeNodesTabFromNewick.py and
#4_1_parseNodes.py.
#Created by David E. Hufnagel on Oct 22, 2012

import sys
from ete2 import Tree

inp = open(sys.argv[1])      #Input newick phylogeny file
out = open(sys.argv[2], "w") #Output parsed newick phylogeny file (subtrees)




def GetDif(bigLst, smallLst):
    difLst = []
    for item in bigLst:
        if item not in smallLst:
            difLst.append(item)
            
    return difLst

#Get species of all descendents
def GetSpecies(tree):
    speSet = set()
    for gene in tree.get_leaves():
        #SITUATION SPECIFIC
        name = gene.name
        spe = SetSpe(name)
        speSet.add(spe)

    return list(speSet)

#Taken directly from 4_1_parseNodes.py
def GetSpe(genes):
    spe = set()
    for gene in genes:
        if gene.startswith("rad") or gene.startswith("R") or gene.startswith("Rr"):
            spe.add("Rr")
        elif (gene.startswith("brs") or gene.startswith("B") or gene.startswith("Br")) and not gene.startswith("Bn"):
            spe.add("Br")
        elif gene.startswith("bns") or gene.startswith("N") or gene.startswith("Bn"):
            spe.add("Bn")
        else:
            print "\nerror here!\n"

    return list(spe)

#Get species of only immediate children
def GetSpeciesChildren(tree):
    speSet = set()
    for gene in tree.get_children():
        #SITUATION SPECIFIC
        name = gene.name
        spe = SetSpe(name)
        speSet.add(spe)

    speLst = list(speSet)
    if "" in speLst:
        speLst.remove("")

    return speLst

def SetSpe(name):
    spe = ""

    if (name.startswith("B") or name.startswith("brs") or name.startswith("Br")) and not name.startswith("Bn"):
        spe = "Br"
    elif name.startswith("R") or name.startswith("rad") or name.startswith("Rr"):
        spe = "Rr"
    elif (name.startswith("N") or name.startswith("bns") or name.startswith("Bn")) and not name.startswith("No"):
        spe = "Bn"
    elif name.startswith("AT"):
        spe = "At"
    elif name.startswith("No"): #just to prevent a printed error
        pass
    else:
        print "ERROR: NEED MORE FUNCTIONALITY"

    return spe

def GetLeftAndRight(tree):
    leftNodes = tree.get_children()[0]
    rightNodes = tree.get_children()[1]

    leftNames = GetNames(leftNodes)
    rightNames = GetNames(rightNodes)

    return leftNames, rightNames

def GetNames(nodes):
    nameLst = []
    for node in nodes:
        nameLst.append(node.name)

    return nameLst


    

#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

treeCnt = 0
for line in inp:
    if not line.startswith("#"):
        treeCnt += 1
        tree = Tree(line.strip("\n"))
        
        #Remove At from tree
        
        ##Remove At from the root node (if it's there)
        topSpe = GetSpeciesChildren(tree) #top level(children of the root nodes) species
        if "At" in topSpe:
            ##Remove the root At
            if tree.get_children()[2].name.startswith("AT"):
                tree.get_children()[2].detach()
                topSpe = GetSpeciesChildren(tree)
                
                ##If there are still "At"s around, get rid of them
                while "At" in topSpe:
                    ##Get rid of top level At leaves
                    if tree.get_children()[0].name.startswith("AT"):
                        tree = tree.get_children()[1]
                        topSpe = GetSpeciesChildren(tree)
                    elif tree.get_children()[1].name.startswith("AT"):
                        tree = tree.get_children()[0]
                        topSpe = GetSpeciesChildren(tree)

                ##Remove internal At leaves (can't use prune function because for some reason it doesn't retain distances the way it should)
                while "At" in GetSpecies(tree):
                    for leaf in tree.get_leaves():
                        if leaf.name.startswith("AT"):
                            parent = leaf.up
                            sis = leaf.get_sisters()[0]
                            if parent == tree:
                                tree = sis
                            else:
                                pD = parent.dist
                                sD = sis.dist
                                leaf.detach()
                                parent.delete()
                                sis.dist = pD + sD
                            break ###Goes back to the while loop instead of the for loop because it's not safe to iterate through a list while deleting members

        #Parse tree
        
        ##Get all non-root nodes
        difLst = GetDif(tree.get_descendants(), tree.get_leaves())
        ##Add the root node
        difLst.append(tree)

        ##Make the "newLineLst" for each node (this is perhaps a poor algorithm, but it is a time saver because it allows me to use the algorithm from 4_1_parseNodes.py directly)
        nodeCnt = 0
        for node in difLst:
            nodeCnt += 1
            
            treeCntStr = str(treeCnt).zfill(4)
            nodeCntStr = str(nodeCnt).zfill(3)
            
            ID = "Tree%s_Node%s" % (treeCntStr, nodeCntStr)
            allSpe = GetSpecies(node) #all species
            leftGenes, rightGenes = GetLeftAndRight(node)
            speStr = ",".join(allSpe)
            leftStr = ",".join(leftGenes)
            rightStr = ",".join(rightGenes)
            newLineLst = [ID, speStr, leftStr, rightStr]

            ##determine whether to keep and therefore to write the line 
            doWrite = True
            speLst = []
            for spe in newLineLst[1].split(","):
                speLst.append(spe)

            if len(speLst) < 3:  ##this test that all species are there AND that there are 3+ genes simultaneously
                doWrite = False
                
            ###test whether any species is represented in both left and right  
            leftGenes = newLineLst[2].split(",")
            rightGenes = newLineLst[3].split(",")
            rightGenes[-1] = rightGenes[-1].strip()

            leftSpe = GetSpe(leftGenes)
            rightSpe = GetSpe(rightGenes)

            testLst = []
            for spe in leftSpe:
                if spe in rightSpe:
                    doWrite = False

            #output info
            if doWrite == True:
                #node.write(format=1, outfile=out)
                out.write(node.write() + "\n")
                out.write("####\n")

                

inp.close()
out.close()
