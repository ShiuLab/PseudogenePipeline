#This script is designed to replace PACID names with gene names in the 1st
#column of a .cds.coord.exon_junc file in the retroPseudogene wrapper
#Created by David E. Hufnagel on April 24, 2013
import sys

inp = open(sys.argv[1])
gff = open(sys.argv[2])
out = open(sys.argv[3], "w")


def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

        
#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through gff and make a dict of key: pacID val: geneName
refDict = {}
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "gene":
            #WARNING: DYNAMIC LINE
            #currGene = "Fves" + ("_".join([lineLst[8].split("Name=")[1].strip("m.g").replace(".","_").split("-")[0], "1", lineLst[8].split("Name=")[1].strip("m.g").replace(".","_").split("-")[1], lineLst[8].split("Name=")[1].strip("m.g").replace(".","_").split("-")[2]])).strip("ene") 
            currGene = lineLst[8].split("Name=")[1].replace(".","_").strip("m.g")
        elif lineLst[2] == "CDS":
            pacID = lineLst[8].split("pacid=")[1]
            if pacID not in refDict:
                refDict[pacID] = currGene

#Go through cds.coord.exon_junc file and replace the PACIDs with the gene names
inp.readline()#Skip the info line
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        pacID = lineLst[0]
        geneName = refDict[pacID]
        out.write("%s\t%s\n" % (geneName, "\t".join(lineLst[1:])))



inp.close()
out.close()
