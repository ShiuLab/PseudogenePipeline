#This script is designed to extract data from an EST file (Rr) and
#GFF files (Br, Al, At) and put it in a karyotype file for use in CIRCOS
#WARNING: HIGHLY SITUATION SPECIFIC
#Created by: David E. Hufnagel on 5-22-2012
#note: must have a "\n" at the end of the original est file
#note: needs to be updated so that the 4th row of a chromosome line should be a raw #
#note: LG determination is a failure
"""Algorithm:
1) Go through Br simplified fasta file and extract size of chromosomes
2) Put that info into chromosome classes
3) Go through Al simplified fasta file and extract size of chromosomes
4) Put that info into chromosome classes
5) Go through At simplified fasta file and extract size of chromosomes
6) Put that info into chromosome classes
7) Go through est and extract the band info for Br, Al and At.
8) Put that info into band classes
9) Estract the EST name, linkage group and cM*1e6 for Rr
10) Put that info into chromosome classes
11) Put that info into band classes
12) Write all info into karyotype file"""

import sys

est = open(sys.argv[1])      #EST file with info about Rr linkage groups and bands and Br, Al, At bands.
RrColor = sys.argv[2]        #the color for Rr chromos
Br = open(sys.argv[3])       #Br simplified fasta file.
BrColor = sys.argv[4]        #the color for Br chromos
Al = open(sys.argv[5])       #Al simplified fasta file.
AlColor = sys.argv[6]        #the color for Al chromos
At = open(sys.argv[7])       #At simplified fasta file.
AtColor = sys.argv[8]        #the color for At chromos
kar = open(sys.argv[9], "w") #karyotype output file.

print "est: ", est
print "RrColor: ", RrColor
print "Br:  ", Br
print "BrColor:  ", BrColor
print "Al:  ", Al
print "AlColor:  ", AlColor
print "At:  ", At
print "AtColor:  ", AtColor
print "kar: ", kar
print





class Chromosome:
    def __init__(self, chromo, coords, color):
        self.chromo = chromo
        if type(coords) == tuple:
            self.coords = coords
        else:
            print "coords should be a tuple"
        self.color = color

    def __str__(self):
        newline = "chr\t-\t%s\t%s\t%s\t%s\t%s\n" % (self.chromo, self.chromo,\
                                                  self.coords[0], self.coords[1],\
                                                  self.color)
        return newline
        
    def UpdateCoords(self, coords):
        if type(coords) == tuple:
            self.coords = coords
        else:
            print "coords should be a tuple"

class Band:
    def __init__(self, chromo, name, coords, color):
        self.chromo = chromo
        self.name = name
        if type(coords) == tuple:
            self.coords = coords
        else:
            print "coords should be a tuple"
        self.color = color

    def __str__(self):
        newline = "band\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.chromo, self.name,\
                                                    self.name, self.coords[0],\
                                                    self.coords[1], self.color)
        return newline

def ExtractChromosomeInfo(fd, color):
    chromosLst = []
    temp = ""
        
    for line in fd:
        coords = (0,0)
        if line.startswith(">"):
            length = 0
            if temp != "":
                chromosLst.append(temp)
            chro = FixName(line[1:-1], fd)
            temp = Chromosome(chro, coords, color)
        #get coords
        else:
            length += len(line[:-1])
            temp.UpdateCoords((0,length))
    #upon exiting the loop add the last chromosome
    else:
        chromosLst.append(temp)
        temp = Chromosome(line[1:-1], coords, color)

    return chromosLst

def ExtractBandInfo(chunk, spe):
    chromo = chunk.split(":")[0]
    coordsStr = chunk.split(":")[1]
    chromo = spe + chromo
    name = "%s:%s" % (chromo, coordsStr)
    coords = (coordsStr.split("-")[0], coordsStr.split("-")[1])
    temp = Band(chromo, name, coords, "black")  #CHANGE COLOR HERE!

    return temp
    
def SaveIntoDict(gene1, gene2, dictX):
    if gene1 not in dictX:
        dictX[gene1] = int(gene2)
    else:
        dictX[gene1] += int(gene2)

def WriteList(chromos):
    for chromo in chromos:
        kar.write(str(chromo))

def FixName(oldName, spe):
    if spe == Br:
        chro = "Br" + oldName
    elif spe == Al:
        chro = "Al" + oldName[-1]
    elif spe == At:
        chro = "At" + oldName
    elif spe == "Rr":
        chro = "Rr" + oldName
    return chro

    



RrRatio = 250000 #a constant for the ratio to multiply the LG by

#get Br, Al, At chromo info
BrChromos = []
BrChromos = ExtractChromosomeInfo(Br, BrColor)

AlChromos = []
AlChromos = ExtractChromosomeInfo(Al, AlColor)

AtChromos = []
AtChromos = ExtractChromosomeInfo(At, AtColor)

#get linkage group coords
est.readline() #skip title line
coordsDict = {}
for line in est:
    lineLst = line.split("\t")
    SaveIntoDict(lineLst[1], lineLst[2], coordsDict)
    
#work with est file (main)
est.seek(0)
est.readline() #skip title line
BrBands = [];AlBands = [];AtBands = []
RrChromos = []; RrBands = []
LGset = set()
for line in est:
    #get Br, Al, At band info
    lineLst = line.split("\t")
    
    BrBand = ExtractBandInfo(lineLst[5], "Br")
    AlBand = ExtractBandInfo(lineLst[6], "Al")
    AtBand = ExtractBandInfo(lineLst[7][:-1], "At")

    BrBands.append(BrBand)
    AlBands.append(AlBand)
    AtBands.append(AtBand)

    #get Rr chromo info
    coords = (0,"1000")  #TO BE DONE MANUALLY (basically I copped out)
    if lineLst[1] not in LGset:  #makes it so only one tuple is kept per LG
        LGset.add(lineLst[1])
        chro = FixName(lineLst[1], "Rr")
        temp = Chromosome(chro, coords, RrColor)
        RrChromos.append(temp)
    
    #get Rr band info
    #Bcoords = (0,0)
    Bcoords = (int(float(lineLst[3]) * RrRatio), int((float(lineLst[3]) * RrRatio) + 1000))
    chro = FixName(lineLst[1], "Rr")
    temp2 = Band(chro, lineLst[0], Bcoords, "black")
    RrBands.append(temp2)

#write all data into kar
WriteList(RrChromos)
WriteList(RrBands)
WriteList(BrChromos)
WriteList(BrBands)
WriteList(AlChromos)
WriteList(AlBands)
WriteList(AtChromos)
WriteList(AtBands)





est.close()
Br.close()
Al.close()
At.close()
kar.close()
