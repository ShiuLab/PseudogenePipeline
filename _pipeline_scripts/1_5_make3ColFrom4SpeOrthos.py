#This script is designed to make a 3 col GO file from a 4 col Orthos file with
#all 4 species (At, Al, Br, Rr)
#Created by David E. Hufnagel on July 3, 2012
"""Algorithm:
1) Go through ref and make a dict of the file
2) Go through inp and extract protein names for the species specified
    3) export new 3 col line by combining line inp line info with related
       refDict info"""

import sys

inp = open(sys.argv[1])      #input 4col orthos file
ref = open(sys.argv[2])      #reference 3col species general GO file
out = open(sys.argv[3], "w") #output species ortho 3col GO file
spe = sys.argv[4]            #2 char species identifier Ex: At  note: 2 is after the spe name for when there are multiple protein names to be concerned with Ex: At2





def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)
        
def OutputData(prot):
    if prot in refDict:
        for annot in refDict[prot]:
            GOID = annot[0]
            GOname = annot[1]
            newLine = "%s\t%s\t%s\n" % (prot, GOID, GOname)
            out.write(newLine)




        
#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through ref and make a dict of the file
refDict = {}    #a dict of key: protName val: GOID, name
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        prot = lineLst[0].split("|")[0]
        SaveIntoDict(prot, (lineLst[1], lineLst[2].strip()), refDict)

#2) Go through inp and extract protein names for the species specified
#    3) export new 3 col line by combining line inp line info with related
#       refDict info
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        
        if spe == "Rr":
            prot = lineLst[3].strip()  #will throw an error if an unnacounted for spe is used because prot will not be defined for use in newLine
        elif spe == "Rr2":
            prots = lineLst[3].strip().split(",")
        elif spe == "At":
            prot = lineLst[0].split(".")[0]

        if spe == "Rr2":
            for prot in prots:
                OutputData(prot)
        else:
            OutputData(prot)




inp.close()
ref.close()
out.close()
