#This script is designed to fragment a genome (At) in a fashion that creates a
#distribution that looks like the Rr genome's distribution
#Created by David E. Hufnagel on July 26, 2012
#Updated on August 2, 2012 to take chunks from the left to right of contigs in Br
#   instead of the middle of the contig (this is what makes it version 2).

"""Algorithm for At:
1) Import Rr contig .size file into a list of sizes
2) Import At fasta file seqs into a list
3) Go through AtSeqLst
    a) Pick a random number from RrSizeLst to define when the loop ends
    b) Run a loop while the length of the sequence is greater than the
       random number just picked
        I) Get the chunk out of At and into a new contig dict with a new name
    c) Output what's left after looping
    d) If the loop was never entered the whole chromosome is the new contig
4) Output sequences with generated names"""

#WARNING: Br algorithm out of date! (as of Aug 2, update)
"""Algorithm for Br:
1) Import Br fasta file into a list of seqs
2) Go through the BrSeqLst
    a) Make a dict of key: size val: [seq,]
    b) Make a list of sizes
3) While Br dict not empty, iterate through
    a) Pick a random number from RrSizeLst to define when the loop ends
    b) Look for a match in Br
        I) If there's a match:
            A) Remove the seq from the val and put the seq in the outLst
            B) If val is now empty, remove the key from the dict
        II) If there's no match:
            A) Look for a bigger contig in Br
                1) If there's a bigger contig:
                    a) Generate a random # between 0 and BrSize-RrSize for the start coord and
                       calculate stop coord
                    b) Pull out seq and add it to outLst
                    c) Make 2 new seqs from what's left and put them in the BrDict
                    d) If the 2 new seqs are longer than 100bp, put them in the BrDict
                    e) Remove the seq from val
                    f) If the val is now empty, remove the key from the Br Dict
                2) If no larger contig is found, count down until you find the closest
                   smaller contig
                    a) Remove the seq from the val and put the seq in the outLst
                    b) If val is now empty, remove the key from the dict
4) Output sequences with generated names"""




import sys, random, bisect

RrSize = open(sys.argv[1])     #Rr genome .size file
genome = open(sys.argv[2])     #input fully intact genome
out = open(sys.argv[3], "w")   #output fragmented genome
spe = sys.argv[4]              #species (At or Br)




def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

#for cases where the randomly generated size was found in the Br dict
#(or where a smaller Br contig was chosen)
def BrRrEqual(size, dictx, outLst, sizeLst):
    #remove the seq from the val and put the seq in the outLst
    seq = dictx[size].pop(0)
    outLst.append(seq)
    
    #if val is now empty, remove the key from the dict and remove the size
    #from the size list
    if dictx[size] == []:
        dictx.pop(size)
        sizeLst.remove(size)

#for cases where the randomly geneated size was not found in the Br dict,
#but a bigger Br contig can be used instead
def BrRrPlus(RrSize, BrSize, dictx, outLst, sizeLst):    
    #pull out seq and add it to outLst
    bigSeq = dictx[BrSize][0]
    newSeq = bigSeq[0:RrSize] #the seq to be the new contig
    outLst.append(newSeq)
    
    #make new seq from what's left and put it in the BrDict
    newSeqR = bigSeq[RrSize:] #the seq created from removing newSeq from bigSeq on the left side

    #if the new seq is longer than 100bp, put it in the BrDict and if it's size
    #isn't already in the size Br size list, add the size to the Br size list
    if len(newSeqR) >= 100:
        SaveIntoDict(len(newSeqR), newSeqR, dictx)
        if not len(newSeqR) in sizeLst:
            bisect.insort(sizeLst, len(newSeqR))

    #remove the seq from val
    dictx[BrSize].pop(0)

    #if val is now empty, remove the key from the dict and remove the size
    #from the size list
    if dictx[BrSize] == []:
        dictx.pop(BrSize)
        sizeLst.remove(BrSize)




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) import Rr contig .size file into a list of sizes
RrSizeLst = []
for line in RrSize:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        RrSizeLst.append(lineLst[1].strip())
            
#for At
if spe == "At":

    #2) import At fasta file seqs into a list
    AtSeqLst = []
    tempStr = ""
    for line in genome:
        if not line.startswith("#"):
            if line.startswith(">"):
                if tempStr != "":
                    AtSeqLst.append(tempStr)
                tempStr = ""
            else:
                tempStr += line.strip()
    #make sure to capture the last sequence on the way out
    else:
        AtSeqLst.append(tempStr)

    #3) go through AtSeqLst...
    outLst = []
    for chromo in AtSeqLst:
        #pick a random number from RrSizeLst to define when the loop ends
        endNum = int(random.choice(RrSizeLst))
        #go through the loop that makes all the new contigs in outLst
        chromoLen = len(chromo)
        isLooped = False
        while len(chromo) > endNum:
            #pick a random number from RrSizeLst to define new contig size
            randSize = int(random.choice(RrSizeLst))
            isLooped = True #simply says that the loop has been entered
            #get the chunk out of At and into a new contig dict with a new name
            newContig = chromo[0:randSize]
            outLst.append(newContig)

            chromo = chromo[randSize:]
        else:
            #output what's left after looping
            if chromo != "":
                if len(chromo) != chromoLen:
                    outLst.append(chromo)

            #if the loop was never entered the whole chromosome is the new contig
            if isLooped == False:
                outLst.append(chromo)


#for Br
elif spe == "Br":
    #Import Br fasta file into a list of seqs
    BrSeqLst = []
    tempStr = ""
    for line in genome:
        if not line.startswith("#"):
            #get the sequence
            if line.startswith(">"):
                if tempStr != "":
                    BrSeqLst.append(tempStr)
                tempStr = ""
            else:
                tempStr += line.strip()
    #make sure to capture the last sequence on the way out
    else:
        BrSeqLst.append(tempStr)

    #Go through the BrSeqLst and make a dict of key: size val: [seq,] and a
    #list of sizes
    BrSeqDict = {}
    BrSizeLst = []
    for seq in BrSeqLst:
        size = len(seq)
        SaveIntoDict(size, seq, BrSeqDict)
        BrSizeLst.append(size)

    #Make BrSizeLst into a sorted list without repeats
    BrSizeLst = list(set(BrSizeLst))
    BrSizeLst.sort()

    #while Br dict not empty, iterate through
    outLst = []
    bigCnt = 1  #just something to be printed to the screen to update progress
    while BrSeqDict != {}:
        #pick a random number from RrSizeLst to define new contig size
        randSize = int(random.choice(RrSizeLst))
        #look for a match in Br
        if randSize in BrSeqDict:
            BrRrEqual(randSize, BrSeqDict, outLst, BrSizeLst)
        else:

            #loop through to work with the closest contig in Br which has a
            #larger size
            biResult = bisect.bisect_right(BrSizeLst, randSize)
            #if no contig is bigger, just take the max value
            if int(biResult) == len(BrSizeLst):
                BrRrPlus(randSize, int(BrSizeLst[biResult - 1]), BrSeqDict, outLst, BrSizeLst)
            #otherwise take the closest bigger contig in Br
            else:
                BrRrPlus(randSize, int(BrSizeLst[biResult]), BrSeqDict, outLst, BrSizeLst)
                isDone = True

        if not bigCnt % 1000:
            print bigCnt
        bigCnt += 1

                
else:
    print "\n***spe error***\n"

        
#4) output sequences with generated names
nameNum = 1
for seq in outLst:
    name = ">%sFake%s\n" % (spe, str(nameNum).zfill(6))
    out.write(name)
    out.write(seq + "\n")
    nameNum += 1




RrSize.close()
genome.close()
out.close()
