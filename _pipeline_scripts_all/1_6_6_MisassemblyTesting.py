#This script is designed to test the rates of putative genomic missasembly.
#More specifically, it takes the distance between the EST match and the
#end of the EST (X), the distance between the matching region on the genome
#and the end of the contig (Y) and the average intron length for the species of interest
#(Z) to calculate whether X >= Y + Z  If true, add a missasemby evidence to
#the contig.  Output the missasembly info for each contig in the format:
#contig   contig   # ESTs   # of evidence     ratio 
#name     size    mapped   for misassembly  ev/ESTs
#RC21015  4639      6            2           0.3333.
#Created by David E. Hufnagel on July 19, 2012
"""Algorithm:
1)  Go through est file and extract coordinates into a dict with a tuple of
    coordinates as the value.
2)  Go through inSize and import file into a dictionary
3)  Go through new est dictionary, get X for each side of the est
    and put the info into a list of ESTs (a new class with name, Xleft,
    and Xright to start with.)
4)  Iterate through ESTs dict
       a) Determine Yleft and Yright
       b) go through ESTList and add Yleft and Yright to the current EST
5)  Calculate if X >= Y + Z for both left and right and add a misassembly
    evidence point to a matching region contig dict where this is true. 
6) Output contig info in the format:
contig   contig   # ESTs   # of evidence     ratio 
 name     size    mapped   for misassembly  ev/ESTs
 RC21015  4639      6            2           0.3333

NOTE, IMPORTANT: Left literally means upstream and Right means downstream in
   all cases in this script
"""


import sys

est = open(sys.argv[1])        #The input unique top-match est file
contigSize = open(sys.argv[2]) #The genomic fasta .size file
estSize = open(sys.argv[3])    #The input est .size file
out = open(sys.argv[4], "w")   #The output file with contig missasembly info
Z = int(sys.argv[5])           #The intron length to be used in the X >= Y + Z calculation





class EST:
    def __init__(self, name, Xleft, Xright, Yleft=-1, Yright=-1):

        #EST name
        if type(name) == str:
            self.name = name
        else:
            print "***Name Type Error!!!***"

        #uptream region length for the EST
        if type(Xleft) == int:
            self.Xleft = Xleft
        else:
            print "***Xleft Type Error!!!***"

        #downtream region length for the EST
        if type(Xright) == int:
            self.Xright = Xright
        else:
            print "***Xright Type Error!!!***"

        #uptream region length for the genomic matching region
        if type(Yleft) == int:
            self.Yleft = Yleft
        else:
            print "***Yleft Type Error!!!***"

        #downtream region length for the genomic matching region
        if type(Yright) == int:
            self.Yright = Yright
        else:
            print "***Yright Type Error!!!***"

    def __str__(self):
        return "name: %s    XYleft: %d,%d    XYright: %d,%d" % \
               (self.name, self.Xleft, self.Yleft, self.Xright, self.Yright)

    #adds a true Yleft and Yright value to the EST class
    #(and would change them if they were already assigned)
    def UpdateY(self, Yleft, Yright):

        #uptream region length for the genomic matching region
        if type(Yleft) == int:
            self.Yleft = Yleft
        else:
            print "***Xleft Type Error!!!***"

        #downtream region length for the genomic matching region
        if type(Yright) == int:
            self.Yright = Yright
        else:
            print "***Xright Type Error!!!***"



class Contig:
    def __init__(self, size, numEST=-1, misEv=-1, ratio=-1.0, ESTs=[]):
        #Contig size
        if type(size) == int:
            self.size = size
        else:
            print "***size Type Error!!!***"

        #The number of ESTs mapped to this contig
        if type(numEST) == int:
            self.numEST = numEST
        else:
            print "***numEST Type Error!!!***"

        #The number of misassembly evidences found (min=0 max=2*numEST)
        if type(misEv) == int:
            self.misEv = misEv
        else:
            print "***misEv Type Error!!!***"

        #The ratio of misEv * 2 /numEST  (we use the *2 to make the max ratio at 1)  
        if type(ratio) == float:
            self.ratio = ratio
        else:
            print "***ratio Type Error!!!***"

        #The ratio of misEv * 2 /numEST  (we use the *2 to make the max ratio at 1)  
        if type(ESTs) == list:
            self.estLst = []
        else:
            print "***ESTs Type Error!!!***"
            
    def __str__(self):
        return "size: %s    numEST: %s    misEv: %s    ratio: %s" % (self.size, self.numEST, self.misEv, self.ratio)
    
    def AddEST(self, anEST):
        self.estLst.append(anEST)
        if self.numEST == -1:
            self.numEST = 1
        else:
            self.numEST += 1

    #Calculates if X >= Y + Z for each EST, and for both left and right and adds
    #a misassembly evidence point to a matching region contig dict where this
    #is true.
    def CalculateMisEvs(self):
        #Left and Right is relative to the EST on its contig
        for eSt in self.estLst:
            #Left
            left = eSt.Xleft >= (eSt.Yleft + Z)
            #Right
            right = eSt.Xright >= (eSt.Yright + Z)

            #add to MisEv
            if left == True:
                if self.misEv == -1:
                    self.misEv = 1
                else:
                    self.misEv += 1
            if right == True:
                if self.misEv == -1:
                    self.misEv = 1
                else:
                    self.misEv += 1

        #calculate and set ratio
        if not (self.misEv == -1 and self.numEST == -1):
            if self.misEv == -1:
                self.ratio = 0
            else:
                self.ratio = float(self.misEv) / (2 * float(self.numEST))
            
        
      

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



print
print "est     :", est
print "contigSize    :", contigSize
print "estSize    :", estSize
print "out      :", out
print "Z        :", Z
print

#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through contig .size file and make a dict of contigs (setting size)
#   key: name val: Contig 
contigDict = {}
for line in contigSize:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        contigDict[lineLst[0]] = Contig(int(lineLst[1].strip()))

estSizeDict = {}
#2) input est .size file into a dict
for line in estSize:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        estSizeDict[lineLst[0].split(" ")[0]] = lineLst[1].strip()

#3) Go through est file, make ESTs in contigs setting all their values
for line in est:
    lineEST = None
    if not line.startswith("#"):
        lineLst = line.split("\t")

        #Get est direction (FR is relative to the est/matching region on its contig)
        frEST = ForwardOrReverse(lineLst[6], lineLst[7]) 
        frMatchRegion =  ForwardOrReverse(lineLst[8], lineLst[9])  

        #Get values from file
        estName = lineLst[0]
        contigName = lineLst[1]
        if frEST == "F":
            Xl = int(lineLst[6]) - 1 #EST.Xleft
            Xr = int(estSizeDict[estName]) - int(lineLst[7])       #EST.Xright
        elif frEST == "R":
            Xl = int(lineLst[7]) - 1 #EST.Xleft
            Xr = int(estSizeDict[estName]) - int(lineLst[6])       #EST.Xright 
        if frMatchRegion == "F":
            Yl = int(lineLst[8]) - 1  #EST.Yleft
            Yr = int(contigDict[contigName].size) - int(lineLst[9]) #EST.Yright
        elif frMatchRegion == "R":
            Yl = int(lineLst[9]) - 1  #EST.Yleft
            Yr = int(contigDict[contigName].size) - int(lineLst[8]) #EST.Yright

        

        ####SPECIFIC TEST###
        if int(lineLst[6]) > int(estSizeDict[estName]) or \
           int(lineLst[7]) > int(estSizeDict[estName]):
            print "EST ERROR AT: %s" % (estName)
        if int(lineLst[8]) > contigDict[contigName].size or \
           int(lineLst[9]) > contigDict[contigName].size:
            print "CONTIG ERROR AT: %s" % (contigName)

        #Put the values into classes
        lineEST = EST(estName, Xl, Xr, Yl, Yr)
        contigDict[contigName].AddEST(lineEST)

        #Calculate X >= Y + Z for each EST and add to misEv adn ratio in Contig
        contigDict[contigName].CalculateMisEvs()

#4) Go through contigs and output info in format:
#contig   contig   # ESTs   # of evidence     ratio 
#name     size    mapped   for misassembly  ev/ESTs
#RC21015  4639      6            2           0.3333.
title = """#contig   contig   # ESTs   # of evidence     ratio 
# name     size    mapped   for misassembly  ev/ESTs\n"""
out.write(title)
for name in contigDict:
    contig = contigDict[name]

    #get num ESTs
    if contig.numEST == -1:
        num = 0
    else:
        num = contig.numEST

    #get num misassembly evidences
    if contig.misEv == -1:
        mis = 0
    else:
        mis = contig.misEv

    #get ratio
    if contig.ratio == -1:
        rat = "-"
    else:
        rat = contig.ratio
        
    newLine = "%s\t%s\t%s\t%s\t%s\n" % (name, contig.size, num, mis, rat)
    out.write(newLine)





est.close()
contigSize.close()
estSize.close() 
out.close()
