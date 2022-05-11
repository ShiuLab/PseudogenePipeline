#This script was designed as a wrapper for getting files for CIRCOS from syntenic
#blocks created by MCScanX
#Created by David E. Hufnagel on Jan 23, 2013
#WARNING: MUST LOAD CIRCOS FIRST

import sys, os
blocks = sys.argv[1]    #input tabular collinearity file
dualSpe = sys.argv[2]   #twoSpecies identifiers Ex: AtBr
genes = sys.argv[3]     #4col file of all genes in the collinearity file (or more)
genomeA = sys.argv[4]   #genome file for the 1st species
genomeB = sys.argv[5]   #genome file for the 2nd species
colors = sys.argv[6]    #the colors to use in the two karyotypes Ex: blue,red
scriptDir = sys.argv[7] #the folder containing all scripts used
scriptDir = "/" + scriptDir.strip("/") #get rid of '/' after scriptDir


if dualSpe == "AtBr" or dualSpe == "BrAt":
    pairs = "AT,Bra"
elif dualSpe == "AtRr" or dualSpe == "RrAt":
    pairs = "AT,Rr"
elif dualSpe == "BrRr" or dualSpe == "BrRr":
    pairs = "Bra,Rr"
else:
    print "ERROR: NEED MORE FUNCTIONALITY"
    sys.exit()
colorLst = colors.split(",")

##prep files
#os.system("python %s/1_2_H_filterCollinearityTab.py ../%s %s %s.%s" % (scriptDir, blocks, pairs, blocks, dualSpe))
os.system("python %s/1_2_H_make4colFromCollinearityTab.py %s.%s %s %sSyntenicBlocks.4colA %sSyntenicBlocks.4colB" % (scriptDir, blocks, dualSpe, genes, dualSpe, dualSpe))
os.system("python %s/1_2_H_get2colFromCollinearityTab.py %s.%s %s %sSyntenicBlocks.2col" % (scriptDir, blocks, dualSpe, genes, dualSpe))
os.system("cat %sSyntenicBlocks.4colA %sSyntenicBlocks.4colB > %sSyntenicBlocks.4colAB" % (dualSpe, dualSpe, dualSpe))
os.system("python %s/circos_links.py %sSyntenicBlocks.4colAB %sSyntenicBlocks.2col %s" % (scriptDir, dualSpe, dualSpe, dualSpe))
os.system("mv %sSyntenicBlocks.2col.links %sSyntenicBlocks.links" % (dualSpe, dualSpe))
os.system("python %s/circos_karyotype.py %s %sSyntenicBlocks.4colA %s %s" % (scriptDir, genomeA, dualSpe, colorLst[0], dualSpe[:2]))
os.system("python %s/circos_karyotype.py %s %sSyntenicBlocks.4colB %s %s" % (scriptDir, genomeB, dualSpe, colorLst[1], dualSpe[-2:]))
os.system("cat %s.karyotype %s.karyotype | grep -v '#' > %sSyntenicBlocks.karyotype" % (genomeA, genomeB, dualSpe))
os.system("mkdir Temp")
os.system("mv %sSyntenicBlocks.4colAB Temp" % (dualSpe))
##run CIRCOS
