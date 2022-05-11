#This script is designed to determine the number of genes in a genome that are
#present in a .collinearity synteny file created by MCScanX
#Created by David E. Hufnagel on April 8, 2013

import sys, os
four = open(sys.argv[1])   #the input .4col genomic genes file
collin = sys.argv[2]       #the input MCScanX .collinearity file




#Go through the 4 column file, make a list of gene names
fourLst = []
for line in four:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        fourLst.append(lineLst[0])

#Tabularize the collinearity file
os.system("python ~/Shiu/Scripts/5_2_2_tabularizeCollinearity.py %s %s.tab" % (collin, collin))

#Go through the tabularized collinearity file and make a set of gene names
tab = open(collin+".tab")
tabLst = []
for line in tab:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
	for pair in lineLst[6:]:
            geneA = pair.split(";")[0]
            geneB = pair.split(";")[1]
	    tabLst.append(geneA)
            tabLst.append(geneB)

#Make sets of the two lists and determine the degree of overlap
fourSet = set(fourLst)
tabSet = set(tabLst)
bothSet = fourSet.intersection(tabSet)
percent = (float(len(bothSet)) / float(len(fourSet))) * 100.0

#Output information
print "%s genes in the genome" % (len(fourSet))
print "%s genes in the collinearity file" % (len(tabSet))
print "%s percent of genes in collinearity file" % (percent)



four.close()
tab.close()
