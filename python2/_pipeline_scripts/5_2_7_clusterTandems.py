#This script is designed to take a tandem gene file created by Gaurav's script,
#find_tandemDuplicates.py in the format: Gene1  Gene2   Distance(# of genes)    Evalue
#and make 3 outputs .clustAll .clustFG .clustPS with tandem clusters (FG means
#functional genes only, PS means pseudogenes only) in the format:
#GroupName  #ofItems  Item1  Item2 Item3...
#Created by David E. Hufnagel on May 16, 2013

import sys
tand = open(sys.argv[1]) #input .tandem file
spe = sys.argv[2]        #species identifier

#output files (see description at top)
outAll = open(sys.argv[1] + ".clustAll", "w")
outFG = open(sys.argv[1] + ".clustFG", "w")
outPS = open(sys.argv[1] + ".clustPS", "w")

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

def RemoveDups(listx):
    newLst = []
    for line in listx:
        if line not in newLst:
            newLst.append(line)

def OutputLst(listx, fd):
    for line in listx:
        fd.write(line)



#Write the users command line prompt on the first line of the output file.
outAll.write("#python %s\n" % (" ".join(sys.argv)))
outFG.write("#python %s\n" % (" ".join(sys.argv)))
outPS.write("#python %s\n" % (" ".join(sys.argv)))

#Go through tandem and make a list of tandem clusters
clustLst = [] # a list of sets where each set is a tandem cluster
cnt = 1
for line in tand:
    if not line.startswith("#"):
        #get features
        lineLst = line.strip().split("\t")
        featA = lineLst[0]
        featB = lineLst[1]
        
        #look for the features in clustLst
        isInClustLst = False #whether one or more features has so far been included in a cluster
        for clust in clustLst:                        
            if featA in clust or featB in clust:
                #print "add both"
                clust.add(featA)
                clust.add(featB)
                isInClustLst = True

        #if the pair should start it's own cluster
        else:    
            if isInClustLst == False:
                newClust = set([featA, featB])
                clustLst.append(newClust)

#Output clustLst with all the proper surrounding info
outAllLst = []
outPSLst = []
outFGLst = []
cnt = 1
for clust in clustLst:
    groupName = "%sTandemClust%s" % (spe, str(cnt).zfill(5))
    numItem = len(clust)

    #output info to outAll and make temporary lists for pseudogenes and genes
    outAllLst.append("%s\t%s" % (groupName, numItem))
    psList = []
    gnList = []
    for item in clust:
        outAllLst.append("\t%s" % (item))
        if item.startswith("Ps"):
            psList.append(item)
        else:
            gnList.append(item)
    outAllLst.append("\n")
    
    #go through gnList and psList and output info to outFG and outPS respectively
    psNumItem = len(psList)
    #If psList is only 0 or 1 items long, don't output it because it's not a cluster
    if psNumItem > 1:
        outPSLst.append("%s\t%s" % (groupName, psNumItem))
        for pseudo in psList:
            outPSLst.append("\t%s" % (pseudo))
        outPSLst.append("\n")

    gnNumItem = len(gnList)
    if gnNumItem > 1:
        outFGLst.append("%s\t%s" % (groupName, gnNumItem))
        for gene in gnList:
            outFGLst.append("\t%s" % (gene))
        outFGLst.append("\n")
    cnt += 1

##remove dups
RemoveDups(outAllLst)
RemoveDups(outPSLst)
RemoveDups(outFGLst)
 
##output output lists
OutputLst(outAllLst, outAll)
OutputLst(outPSLst, outPS)
OutputLst(outFGLst, outFG)




                
tand.close()
outAll.close()
outFG.close()
outPS.close()
