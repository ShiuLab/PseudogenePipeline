#This script is designed to fragment a genome (At) in a fashion that creates a
#distribution that looks like the Rr genome's distribution
#Created by David E. Hufnagel on July 26, 2012
#Updated on August 2, 2012 to take chunks from the left to right of contigs in Br
#   instead of the middle of the contig (this is what makes it version 2).
#Updated on August 7, 2012 to make an output file with information about how
#   the new contigs relate to the old ones. (this is what makes it version 3.)

#WARNING: At algorithm out of date! (as of Aug 7 update)
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
ref = open(sys.argv[4], "w")   #output reference file with information about how the new contigs relate to the old ones




def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)

def ImportFasta(fasta):  #fasta is the fasta file name
    fastaDict = {}
    currName = ""   #the current name associated with a seq
    seq = ""        #the current seq to be built up for each seq line
    for line in fasta:
        if not line.startswith("#"):
            if line.startswith(">"):
                if currName != "":
                    fastaDict[currName] = seq
                seq = ""
                currName = line.strip().strip(">")
            else:
                seq += line.strip()
                
    #get the last seq on the way out
    else:
        fastaDict[currName] = seq

    return fastaDict




#Write the users command line prompt on the first line of the output files.
out.write("#python %s\n" % (" ".join(sys.argv)))
ref.write("#python %s\n" % (" ".join(sys.argv)))

#1) import Rr contig .size file into a list of sizes
RrSizeLst = []
for line in RrSize:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        RrSizeLst.append(lineLst[1].strip())
            
#2) import At fasta file seqs into a list
fastaDict = ImportFasta(genome)

#3) go through AtSeqLst...
outLst = []  #format [(seq, oldChName, oldChStartCoord, oldChStopCoord),]
for name, seq in fastaDict.items():
    #pick a random number from RrSizeLst to define when the loop ends
    endNum = int(random.choice(RrSizeLst))
    #go through the loop that makes all the new contigs in outLst
    chromoLen = len(seq)
    isLooped = False
    startOld = 0
    while len(seq) > endNum:
        #pick a random number from RrSizeLst to define new contig size
        randSize = int(random.choice(RrSizeLst))
        isLooped = True #simply says that the loop has been entered
        #get the chunk out of At and into a new contig dict with a new name
        newContig = seq[0:randSize]
        outLst.append((newContig, name, startOld + 1, startOld+randSize))

        seq = seq[randSize:]
        startOld += randSize
    else:
        #output what's left after looping
        if seq != "":
            if len(seq) != chromoLen:
                outLst.append(seq)

        #if the loop was never entered the whole chromosome is the new contig
        if isLooped == False:
            outLst.append(seq)

#4) output sequences with generated names and output info into ref file
nameNum = 1
for tup in outLst:
    seq = tup[0]
    name = ">AtFake%s\n" % (str(nameNum).zfill(6))
    out.write(name)
    out.write(seq + "\n")
    nameNum += 1

    refLine = "%s:%s-%s\t%s:%s-%s\n" % (name.strip().strip(">"), 1, len(seq), tup[1], tup[2], tup[3])
    ref.write(refLine)




RrSize.close()
genome.close()
out.close()
ref.close()
