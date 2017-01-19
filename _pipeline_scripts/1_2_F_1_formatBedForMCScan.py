#This script is designed to take a .bed file for all the genes in Br and change
#The chromosome names so they conform to MCScanX input requirements  AKA:
#BrA10 -> Br000010 and BrScaffold000096 -> Br000096
#Created by David E. Hufnagel on Nov 1, 2012
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])      #input .bed file with normal names
out = open(sys.argv[2], "w") #output .bed file with reformed names
ref = open(sys.argv[3], "w") #output tabular reference file with 1st col: newName 2nd col: oldName




#Write the users command line prompt on the first line of the output files.
out.write("#python %s\n" % (" ".join(sys.argv)))
ref.write("#python %s\n" % (" ".join(sys.argv)))

refLst = [] #keeps track of what names are already written into the ref file
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip("\n").split("\t")

        #For the At chromosomes
        if "AtC" in lineLst[0]:
            newName = "".join(lineLst[0].split("Chr"))
        #For the Br chromosomes
        elif "BrA" in lineLst[0]:
            newName = lineLst[0].split("A")[0] + lineLst[0].split("A")[1]#.strip("0").zfill(6)
        #For the Br scaffolds
        elif "BrS" in lineLst[0]:
            newName = "".join(lineLst[0].split("Scaffold"))
        #For the Rr scaffolds
        elif "RrC" in lineLst[0]:
            newName = "".join(lineLst[0].split("C"))
        #For the Al scaffolds
        elif "scaffold_" in lineLst[0]:
            newName = "Al" + lineLst[0].split("scaffold_")[1]

        #output info into out
        outLine = "%s\t%s\t%s\t%s\n" % (newName, lineLst[1], lineLst[2], lineLst[3])
        out.write(outLine)

        #output unique info into ref
        if newName not in refLst:
            refLst.append(newName)
            refLine = "%s\t%s\n" % (newName, lineLst[0])
            ref.write(refLine)
            




inp.close()
out.close()
ref.close()
