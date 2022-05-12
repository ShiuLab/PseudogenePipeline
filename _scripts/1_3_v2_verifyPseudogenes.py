#This script is designed to verify pseudogenes based on their distance
#to the end of contigs and the length of the overhang of its corresponding
#paralagous functional gene.  More specifically, when pseudogenes do not contain
#stop codons or frameshits it takes the distance between the pseudogene and the
#end of the contig (X), the distance between the matching region on the protein
#and the end of the protein (Y) and the 95th percentile intron length for the
#species being tested (Z) to calculate whether X >= Y + Z  If the pseudogene
#satisfies this criteria (X is greater), it is output into a high confidence
#file of each format (.disable_count, .real4col, .fa).
#Created by David E. Hufnagel on June 19, 2012
"""Algorithm:
1)  Go through inFour and extract pseudogene coordinates into a dict with a
    tuple of coordinates as the value.
2)  Make a new Dict of pseudogenes that do NOT have any stop codons or frameshifts.
3)  Go through inSize and import file into a dictionary
4)  Go through new pseudogene dictionary, get X for each side of the pseudogene
    and put the info into a list of Pseudogenes (a new class with name, Xleft,
    and Xright to start with.)
5)  Go through inProtSize and import names and sizes of proteins into a dictionary
6)  Iterate through questionable pseudogene names dict
       a) Determine Yleft and Yright
       b) go through pseuList and add Yleft and Yright to the current pseudogene
7)  Calculate if X >= Y + Z for both left and right and make a badLst of the
    names that didn't make the cut
8)  Go through inDis and filter to make outDis
9)  Go through inFour and filter to make outFour
10) Go through inFasta and filter to make outFasta

NOTE, IMPORTANT: Left literally means upstream and Right means downstream in
   all cases in this script
"""

import sys

inDis = open(sys.argv[1])         #The disable_count format pseudogene input file for filtering       (.fullyFiltered.disable_count.RMfilt)
inFour = open(sys.argv[2])        #The four column pseudogene input file for filtering                (.real4col.RMfilt)
inFasta = open(sys.argv[3])       #The fasta format pseudogene input file for filtering               (.real4col.mod.fa.RMfilt)
inSize = open(sys.argv[4])        #The genomic fasta .size file                                       ([genome].fa.size)                      
inProtSize = open(sys.argv[5])    #The protein fasta .size file                   		      ([protein].fa.size)
outDis = open(sys.argv[6], "w")   #The output high confidence pseudogenes in the disable_count format (.fullyFiltered.disable_count.RMfilt.hiConf)
outFour = open(sys.argv[7], "w")  #The output high confidence pseudogenes in the four column format   (.real4col.RMfilt.hiConf)
outFasta = open(sys.argv[8], "w") #The output high confidence pseudogenes in the fasta format         (.real4col.mod.fa.RMfilt.hiConf)
Z = int(sys.argv[9])              #The intron length to be used in the X >= Y + Z calculation          (intron length)

class Pseudogene:
    def __init__(self, name, Xleft, Xright, Z, Yleft=0, Yright=0):

        #pseudogene name
        if type(name) == str:
            self.name = name
        else:
            print("***Name Type Error!!!***")

        #uptream region length for the pseudogene
        if type(Xleft) == int:
            self.Xleft = Xleft
        else:
            print("***Xleft Type Error!!!***")

        #downtream region length for the pseudogene
        if type(Xright) == int:
            self.Xright = Xright
        else:
            print("***Xright Type Error!!!***")

        #estimated intron length
        if type(Z) == int:
            self.Z = Z
        else:
            print("***Z Type Error!!!***")

        #uptream region length for the protein matching region
        if type(Yleft) == int:
            self.Yleft = Yleft
        else:
            print("***Yleft Type Error!!!***")

        #downtream region length for the protein matching region
        if type(Yright) == int:
            self.Yright = Yright
        else:
            print("***Yright Type Error!!!***")

    def __str__(self):
        return "name: %s    XYleft: %d,%d    XYright: %d,%d    Z: %d" % \
               (self.name, self.Xleft, self.Yleft, self.Xright, self.Yright, self.Z)

    #adds a true Yleft and Yright value to the pseudogene class
    #(and would change them if they were already assigned)
    def UpdateY(self, Yleft, Yright):

        #uptream region length for the protein matching region
        if type(Yleft) == int:
            self.Yleft = Yleft
        else:
            print("***Xleft Type Error!!!***")

        #downtream region length for the protein matching region
        if type(Yright) == int:
            self.Yright = Yright
        else:
            print("***Xright Type Error!!!***")

    #Calculates if X > Y + Z and returns the boolean result
    def IsHighConf(self):
        #Left and Right is relative to the pseudogene on its contig
        isLeft = (self.Xleft >= (self.Yleft + Z))
        isRight = (self.Xright >= (self.Yright + Z))
        isBoth = (isLeft and isRight)

        return isBoth

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

def ForwardOrReverse(start, stop):
    if int(start) < int(stop):
        return "F"
    elif int(start) > int(stop):
        return "R"
    else:
        return "***FR error!!!***"

"""
print()
print("inDis      :", inDis)
print("inFour     :", inFour)
print("inFasta    :", inFasta)
print("inSize     :", inSize)
print("inProtSize :", inProtSize)
print("outDis     :", outDis)
print("outFour    :", outFour)
print("outFasta   :", outFasta)
print("Z          :", Z)
print()

#writes the users command line prompt on the first line of the output file.
outDis.write("#python %s\n" % (" ".join(sys.argv)))
#writes the users command line prompt on the first line of the output file.
outFour.write("#python %s\n" % (" ".join(sys.argv)))
#writes the users command line prompt on the first line of the output file.
outFasta.write("#python %s\n" % (" ".join(sys.argv)))
"""

#1) Go through inFour and extract pseudogene coordinates into a dict with a
#   tuple of coordinates as the value.
pseuDict = {}
for line in inFour:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[0], (lineLst[2], lineLst[3][:-1]), pseuDict)

quesPseuDict = {} #questionable pseudogene dictionary to be sent through the filter
#2) Make a new Dict of pseudogenes that do NOT have any stop codons or frameshifts.
for key in pseuDict:
    val = pseuDict[key]
    code = key.split(";")[-1].split("|")

    #questionable pseudogenes
    if code[0] == "0" and code[1] == "0" and code[2] == "0" and code[3] == "0":
        SaveIntoDict(key, val[0], quesPseuDict) #val[0] is because its always a list of 1

#3) Go through inSize and import file into a dictionary
sizeDict = {}
for line in inSize:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[0], lineLst[1][:-1], sizeDict)

#4) Go through new pseudogene dictionary, get X for each side of the pseudogene
#   and put the info into a list of Pseudogenes (a new class with name, Xleft,
#   and Xright to start with.)
pseuList = []
for key in quesPseuDict:
    val = quesPseuDict[key]
    #fr is relative to the pseudogene on its contig
    fr = ForwardOrReverse(val[0][0], val[0][1]) #starts with 0 because it's always a 1 member list because of the way SaveIntoDict works
    chro = "|".join(key.split(";")[1].split("|")[:-1])

    if fr == "F":
        Xl = int(val[0][0]) - 1
        Xr = int(sizeDict[chro][0]) - int(val[0][1]) #starts with 0 because it's always a 1 member list because of the way SaveIntoDict works
    elif fr == "R":
        Xl = int(sizeDict[chro][0]) - int(val[0][0]) #starts with 0 because it's always a 1 member list because of the way SaveIntoDict works
        Xr = int(val[0][1]) - 1 
    else:
        print("***FR error 2***")
        
    newPseu = Pseudogene(key, Xl, Xr, Z)
    pseuList.append(newPseu)

# Modified by Nick Panchy
#5a) Create a dictioanry of protein sequence sizes
protein_size = {}
for line in inProtSize:
    if not line.startswith("#"):
        split_line = line.strip().split("\t")
        protein_size[split_line[0]] = int(split_line[1])

#6a)Iterate through questintable pseudogenes
cnt = 0
for key in quesPseuDict:
    val = quesPseuDict[key]
    prot_name = key.split(";")[0]#.split("|")[0]
    pseuProt = tuple(key.split(";")[2].split(":")[0].split("-")) #pseudogene coords on the protein
    size = protein_size[prot_name]
    Yl = (int(pseuProt[0]) - 1) * 3
    Yr = (size - int(pseuProt[1]) -1)* 3
    # go through pseuList and add Yleft and Yright to the current pseudogene
    for pseu in pseuList:
        if pseu.name == key:
            pseu.UpdateY(Yl, Yr)

badLst = []
#7) Calculate if X >= Y + Z for both left and right and make a badLst of the
#   names that didn't make the cut
for pseu in pseuList:
    isHC = pseu.IsHighConf()
    if isHC == False:
        badLst.append(pseu.name)
    #print pseu

#8) Go through inDis and filter to make outDis
write = True
for line in inDis:
    if line.startswith("#"):
        lineLst = line.split(" ")
        name = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:])[:-1])[1:]
        if not name in badLst:
            write = True
            outDis.write(line)
        else:
            write = False

    else:
        if write == True:
            outDis.write(line)

#9) Go through inFour and filter to make outFour
inFour.seek(0)
for line in inFour:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if not lineLst[0] in badLst:
            outFour.write(line)

#10) Go through inFasta and filter to make outFasta
write = True
for line in inFasta:
    if line.startswith(">"):
        if not line[1:-1] in badLst:
            write = True
            outFasta.write(line)
        else:
            write = False
    else:
        if not line.startswith("#"):
            if write == True:
                outFasta.write(line)

inDis.close()
inFour.close()
inFasta.close()
inSize.close()
inProtSize.close()
outDis.close()
outFour.close()
outFasta.close()
