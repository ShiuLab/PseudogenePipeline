#This script is designed to make a tabular file from a Newick phylogeny file.
#The output file is all nodes in the trees in the format:
#ID  AllSpecies LeftGenes  RightGenes

import sys
from ete2 import Tree

inp = open(sys.argv[1])      #Input newick phylogeny file
out = open(sys.argv[2], "w") #Output tabular file




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

    if name.startswith("B") or name.startswith("brs") or name.startswith("Br"):
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
    if not line.startswith("#p"):
        if not line.startswith("####"):
            treeCnt += 1
            tree = Tree(line.strip("\n"))
            #remove At from the root node (if it's there)
            topSpe = GetSpeciesChildren(tree) #top level(children of the root nodes) species
            if "At" in topSpe:
                #remove the root At
                if tree.get_children()[2].name.startswith("AT"):
                    tree.get_children()[2].detach()
                    topSpe = GetSpeciesChildren(tree)
                    
                    #if there are still "At"s around, get rid of them
                    while "At" in topSpe:
                        #Get rid of top level At leaves
                        if tree.get_children()[0].name.startswith("AT"):
                            tree = tree.get_children()[1]
                            topSpe = GetSpeciesChildren(tree)
                        elif tree.get_children()[1].name.startswith("AT"):
                            tree = tree.get_children()[0]
                            topSpe = GetSpeciesChildren(tree)

                    #Remove internal At leaves (can't use prune function because for some reason it doesn't retain distances the way it should)
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
                                break #goes back to the while loop instead of the for loop because it's not safe to iterate through a list while deleting members

            #get all non-root nodes
            difLst = GetDif(tree.get_descendants(), tree.get_leaves())
            #add the root node
            difLst.append(tree)
            #output info
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
                
                newLine = "%s\t%s\t%s\t%s\n" % (ID, speStr, leftStr, rightStr)
                out.write(newLine)
            



inp.close()
out.close()
