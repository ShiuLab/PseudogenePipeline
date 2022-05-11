#This script is designed to put all phylogenetic groups (newick format) into
#cases for analysis, and therefore counts the cases and the amount of things not
#in cases.
#WARNING: This is specific to the B. nigra B. rapa hybrid problem.
#Created by David E. Hufnagel on Sep 24, 2012

import sys
from ete2 import Tree

inp = open(sys.argv[1]) #The input newick format phylogeny file
#results are printed to the screen




#Removes the numbers from the phylogeny so only gene names remain
def RemoveNumbers(string):
    keep = True
    newStr = ""
    for char in string:
        if char == ":":
            keep = False
        elif char == "," or char == ")":
            keep = True
            newStr += char
        elif keep == True:
            newStr += char

    return newStr

#Makes the gene names all strings when hit by eval.  Always immediately follows RemoveNumvers()
def MakeNamesStrings(string):
    inStr = False
    newStr = ""
    for char in string:
        #start string
        if char not in "()," and inStr == False:
            inStr = True
            newStr += "'"
        #end string
        elif char in "()," and inStr == True:
            inStr = False
            newStr += "'"
        newStr += char

    return newStr

#get simple list of genes in tree. Cnt should start at 0
def GetListGenes(string):
    inStr = False
    tempStr = ""
    geneLst = []
    for char in string:
        if inStr == True:
            #end string
            if char in "(),":
                inStr = False
                geneLst.append(tempStr)
                tempStr = ""
            #continuing string
            else:
                tempStr += char
        else:
            #start string
            if char not in "(),":
                inStr = True
                tempStr += char

    return geneLst

#Reduces the tree by combining genes of the same species when they are the closest in the tree
def ReduceTreeSimple(tree, treeStr):
    treeLst = GetListGenes(treeStr)  #I didn't just pass in the list because this is a temporary list

    #do a very simple check: if two of the same species aren't next to each other in the treeLst, they cannot be reduced.
    doProcess = False #whether to go past this test to the full analysis
    lastStart = ""
    for gene in treeLst:
        #HIGHLY SITUATION SPECIFIC
        if gene[:3] == lastStart:
            doProcess = True
        lastStart = gene[:3]

    #do actual reduction        

#using a tree made from the tree library
def ReduceTreeSimple2(root):
    print root
    stillMerging = True
    while stillMerging == True:
        print
        toMerge = []
        stillMerging = False #will be reset to true if necessary
        #find the leaf to merge
        for node in root.iter_leaves():
            print node
            print node.get_sisters()
            sisName = node.get_sisters()[0].name
            if sisName != "NoName": #it's a single member list
                if node.name[:3] == sisName[:3]:
                    if node.name not in toMerge and sisName not in toMerge:
                        toMerge.append(node.name)
                        stillMerging = True
    ##        print node.name
    ##        print "sis: ", node.get_sisters()
    ##        print
    ##
        #merge it
        for nodeN in toMerge:
            #get the dist of the common ancestor
            node1 = root.search_nodes(name = nodeN)[0]
            node2 = node1.get_sisters()[0]
            ancestorD = root.get_common_ancestor(nodeN, node2.name).dist

            #set the dist to one node and delete the other.  Also rename the remaining node
            #print node1.name
            #print node1.dist
            node1.name = nodeN + ";" + node2.name
            node1.dist = ancestorD
            node2.detach()
            #print node1.name
            #print node1.dist

        print root
    sys.exit()
    return root

#Break up tree into units of 3 when possible
def ReduceTreeByThree(root):
    #print root

    #determine if the tree is reducable
    isReducable = False
    if not len(root.get_leaves()) % 3 and len(root.get_leaves()) != 3:
        #for 6
        if len(root.get_leaves()) == 6:
            #print 6
            isReducable = True
            for child in root.get_children():
                if len(child.get_leaves()) != 3:
                    isReducable = False

        #for 9
        elif len(root.get_leaves()) == 9:
            #print 9
            #see if the children of the root have 6 and 9 members
            isThree = False
            isSix = False
            for child in root.get_children():
                if len(child.get_leaves()) == 3:
                    isThree = True
                elif len(child.get_leaves()) == 6:
                    isSix = True
            if isThree == True and isSix == True:
                #see if the children of the 6-member child each have 3 members
                isReducable = True
                for child in root.get_children():
                    if len(child.get_leaves()) == 6:
                        for baby in child.get_children():
                            if len(baby.get_leaves()) != 3:
                                isReducable = False
            
        else:
            print "need more capabilities"

    #reduce it
    if isReducable == True:
        pass
        #make the new trees
        #delete the old ones
        
    
    
    #sys.exit()
        



temp = 0
treeLst = []
for line in inp:
    if not line.startswith("#p"):
        if line.startswith("####"):
            pass
        else:
            tree = Tree(line.strip("\n"))
            temp += 1
##            if temp == 11:
##                tree = ReduceTreeSimple2(tree)
##            if temp == 31:
##                print line
            #print line
            #if temp == 100:
            print tree
            ReduceTreeByThree(tree)
            treeLst.append(tree)
            
            newLineTemp = RemoveNumbers(line.strip("\n").strip(";"))
            geneLst = GetListGenes(newLineTemp)
            newLine = MakeNamesStrings(newLineTemp)
            newLst = eval(newLine)

            ### test cases ###
            #Case 1: 3 genes, one of each species (there are none of these
            isOneOfEach = True
            tempLst = []
            case1Cnt = 0
            for gene in geneLst:
                #HIGHLY SITUATION SPECIFIC
                if gene[:3] in tempLst:
                    isOneOfEach = False
                    break
                else:
                    tempLst.append(gene[:3])

            if len(geneLst) == 3 and isOneOfEach == True:
                case1Cnt += 1

            #Case 2: reducible to case 1
            ReduceTreeSimple(newLst, newLineTemp)
            
    
