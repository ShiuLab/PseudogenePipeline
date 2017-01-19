#This script is designed to make a links file for CIRCOS from an MCscan .aligns
#output file.
#Created by David E. Hufnagel on 11-11-11
#note: is capable of dealing with coordinates in + or - direction
#Algorithm:
"""
1. Skip Header lines
2. Get reference pairs (should be one for each block, in the middle of the block)
3. import .bed file into a list of lists (so I can iterate through it more than once)
4. Get coordinate data for those pairs (from the .bed file) in the same place as the pairs,
     in a Reference pair class.
5. Make a unique value for each pair and write the reference pair info in the
   CIRCOS links format in the ouput file.
"""





class RefPair:

    #                     int       str         str       int          int         int        int          str       str
    def __init__( self, num = 0, name1 = "", name2 = "", coor1a = 0, coor1b = 0, coor2a = 0, coor2b = 0, chr1 = "", chr2 = ""):
        self.mNum = num        #The number representing the Reference gene pair
        self.mName1 = name1    #The name of the first gene in the gene pair
        self.mName2 = name2    #The name of the second gene in the gene pair
        self.mCoor1a = coor1a  #The first coordinate for the first gene
        self.mCoor1b = coor1b  #The second coordinate for the first gene
        self.mCoor2a = coor2a  #The first coordinate for the second gene
        self.mCoor2b = coor2b  #The second coordinate for the second gene
        self.mChr1 = chr1      #The name of the chromosome the first gene is located on
        self.mChr2 = chr2      #The name of the chromosome the second gene is located on

    def __str__( self):
        return "%d:: %s:%d, %d %s:%d, %d" % \
               (self.mNum, self.mName1, self.mCoor1a, self.mCoor1b,\
                self.mName2, self.mCoor2a, self.mCoor2b)



class Block:

    #                   int      int                int       int                 int      list[RefPairs]
    def __init__( self, num = 0, min1 = 1000000000, max1 = 0, min2 = 1000000000, max2 = 0, pairLst = []):
        self.mNum = num               #The block number (index)
        self.mMin1 = min1             #The min value for the first set of genes
        self.mMax1 = max1             #The max value for the first set of genes
        self.mMin2 = min2             #The min value for the second set of genes
        self.mMax2 = max2             #The max value for the second set of genes
        self.mPairLst = pairLst       #The list of reference pairs this block contains
        self.mNumPairs = len(pairLst) #The number of pairs in the pair list
        #mChromo1 and 2 are based on the tested assumption that all genes in a block will be from the same chromosome
        self.mChromo1 = self.mPairLst[0].mChr1 #The chromosome for the first set of genes
        self.mChromo2 = self.mPairLst[0].mChr2 #The chromosome for the second set of genes
        
    def __str__( self):
        return "%d:: min1: %d, max1: %d, min2: %d, max2: %d, NumPairs: %d" % \
               (self.mNum, self.mMin1, self.mMax1, self.mMin2, self.mMax2, self.mNumPairs)
            
 

def SetCoordsAndChromos(pair): #pair is of type RefPair

    for lineLst in bedLst:
        #set coords for name1
        if lineLst[3].strip() == pair.mName1:
            if int(lineLst[1]) < int(lineLst[2]): # if + directionality
                pair.mCoor1a = int(lineLst[1])
                pair.mCoor1b = int(lineLst[2])
            elif int(lineLst[1]) > int(lineLst[2]): # if - directionality
                pair.mCoor1a = int(lineLst[2])
                pair.mCoor1b = int(lineLst[1])
            else:
                print "\nERROR: UNEXPECTED COORDINATES\n"

            #set chromosome name1
            pair.mChr1 = lineLst[0]

        #set coords for name2
        elif lineLst[3].strip() == pair.mName2:        
            if int(lineLst[1]) < int(lineLst[2]): # if + directionality
                pair.mCoor2a = int(lineLst[1])
                pair.mCoor2b = int(lineLst[2])
            elif int(lineLst[1]) > int(lineLst[2]): # if - directionality
                pair.mCoor2a = int(lineLst[2])
                pair.mCoor2b = int(lineLst[1])
            else:
                print "\nERROR: UNEXPECTED COORDINATES\n"

            #set chromosome name2
            pair.mChr2 = lineLst[0]

def ImportBed():

    newLst = []
    for line in bed:
        lineLst = line.split("\t")
        newLst.append(lineLst)

    return newLst

def MakeBlocks(tempLst): # pass in a list of lists of reference pairs
    blockLst = [] #a list of blocks of type Block
    cnt = 0
    for block in tempLst:
        coorList1 = []
        coorList2 = []
        for pair in block:
            coorList1.append(pair.mCoor1a)
            coorList1.append(pair.mCoor1b)
            coorList2.append(pair.mCoor2a)
            coorList2.append(pair.mCoor2b)

        min01 = min(coorList1)
        max01 = max(coorList1)
        min02 = min(coorList2)
        max02 = max(coorList2)

        tempBlock = Block(num = cnt, min1 = min01, max1 = max01,\
                          min2 = min02, max2 = max02, pairLst = block)
        blockLst.append(tempBlock)
        
        cnt += 1
    return blockLst

def AddZeroes(num): # 3  -> "003"
    tempNum = str(num)
    while len(tempNum)<3:
        tempNum = "0" + tempNum    
    
    newNum = tempNum
    return newNum

import sys

aligns = open(sys.argv[1])   #the name of the .aligns input file (for syntenic blocks)
bed = open(sys.argv[2])      #the name of the .bed input file (for coordinate data)
out = open(sys.argv[3], "w") #the name of the output file
outTest = open("testAssumption.txt", "w") #the output file for testing that each block contains only one kind of chromosome/scaffold for each species

#writes the users command line prompt on the first line of the output file.
out.write('#python %s\n'%(' '.join(sys.argv)))



    

## add gene pairs ##
refCnt = 1      #The number representing the current number of reference pairs
thisBlock = []  #The temporary list of RefPairs in one block
tempLst = []    #big temporary list of "thisBlock" lists (of RefPairs)
isHeader = True #boolean to state whether we are currently iterating in the informational header
for line in aligns:
    if not isHeader: #not in the header
        if line.startswith("## A"): #a block tile line
            refCnt += 1
            tempLst.append(thisBlock)
            thisBlock = []
        else: #actual pairwise data
            lineLst = line.split("\t")
            #if refCnt == 0:
                #print lineLst
            #set the reference pair object
            temp = RefPair(num = refCnt, name1 = lineLst[1], name2 = lineLst[2])

            # Get coordinates
            #SetCoords(temp)
            
            thisBlock.append(temp)
            
    else: 
        if line.startswith("## A"): #The first line of the code (not header)
            isHeader = False



## import .bed file into a list of lists ##
bedLst = ImportBed()
 
## add coordinates and chromosomes from .bed file##
cnt = 0
for block in tempLst:
    for pair in block:
        cnt += 1
        #if cnt < 1000:
        SetCoordsAndChromos(pair)
        if not (cnt % 1000):
            print cnt
        #print pair
        """print pair.mChr1, ":  ", pair
        print pair.mChr2, ":  ", pair
        print"""

#makes sure every block is part of just one chromosome and one scaffold
cnt = 0
for block in tempLst:
    tmp1 = ""
    tmp2 = ""
    flag1 = False
    flag2 = False
    for pair in block:
        toWrite = "%d::  %s : %s\n" % (cnt, pair.mChr1, pair.mChr2)
        outTest.write(toWrite)
        """if tmp1 != pair.mChr1:
            if flag1 == True:
                print "complain 1"
            flag1 = True
        if tmp2 != pair.mChr2:
            if flag2 == True:
                print "complain 2"
            flag2 = True
        
        tmp1 = pair.mChr1
        tmp2 = pair.mChr2"""
    cnt += 1
        
        

## Get Min and Max Coords for a block and make blockLst##
blockLst = MakeBlocks(tempLst)

## output data into links file ##
cnt = 0
for block in blockLst:
    cnt2 = AddZeroes(cnt)
    name = "AthAly" + cnt2

    outStr1 = "%s\t%s\t%d\t%d\n" % (name, block.mChromo1, block.mMin1, block.mMax1)
    outStr2 = "%s\t%s\t%d\t%d\n" % (name, block.mChromo2, block.mMin2, block.mMax2) 
    out.write(outStr1)
    out.write(outStr2)
    
    cnt += 1
    


#close file pipelines
aligns.close()
bed.close()
out.close()
outTest.close()
