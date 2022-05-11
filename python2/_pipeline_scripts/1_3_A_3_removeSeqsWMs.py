#This script is designed to go through a m8 blast file and remove all
#seqs with "M"s from masking
#Created by David E. Hufnagel on Aug 6, 2012

import sys
blast = open(sys.argv[1])    #input m8 blast file
fasta = open(sys.argv[2])    #fasta file containing seqs with "M"s
out = open(sys.argv[3], "w") #output m8 blast file without seqs with "M"s



def OrderCoords(coor1, coor2):
    if type(coor1) != int or type(coor2) != int:
        print "***TYPE ERROR***"
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        print "***ERROR HERE***"

    return small, big

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



#go through fasta and inport the file into a dict of key: name val: seq
fastaDict = ImportFasta(fasta)

#go through blast and output seqs without M's in them
for line in blast:
    if not line.startswith("#"):
        lineLst = line.split("\t")

        #get seq from coords
        chromo = lineLst[1]
        startInd, endInd = OrderCoords(int(lineLst[8]) - 1, int(lineLst[9]) - 1)
        seq = fastaDict[chromo][startInd:endInd+1]

        #if seq is without "M"s write the blast line into the output
        if "M" not in seq:
            out.write(line)
 



blast.close()
fasta.close()
out.close()
