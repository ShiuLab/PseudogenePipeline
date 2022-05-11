#This script is designed to convert a pseudogene file
#(.4col.true.condensed.noCM),  into nucleotide fasta format
"""Algorithm:
1) Go through pseudogene file and make a dict of names and coordinates
2) Go through genome file:
    a) extract genome seq related to coordinates
    b) make dicts of names and seqs from extraction
3) Output dict into new file"""
#WARNING: NEVER COMPLETED

import sys

pseu = open(sys.argv[1])     #input pseudogene file (.4col.true.condensed.noCM)
genome = open(sys.argv[2])   #input genome file
out = open(sys.argv[3], "w") #output file (.fullyFiltered.fa)



#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through pseudogene file and make a dict of names and coordinates
coordDict = {}
pseu.readline()
for line in pseu:
    lineLst = line.split("\t")
    coords = ((lineLst[2], lineLst[3]))
    chromo = lineLst[0].split(";")[1].split("|")[0]
    if lineLst[0] not in coordDict:
        coordDict[lineLst[0]] = [coords, chromo]  #the value is a list with a tuple (the coords) and the chromo string in that order Ex: [(456639, 456789), Chr03]
    else:
        print "problem!!!"

#2) Go through genome file:
#    a) extract genome seq related to coordinates
#    b) make dicts of names and seqs from extraction
read = False
tempLen = 0
lastTempLen = 0
finalDict = {}
for key in coordDict:
    dictChromo = coordDict[key][1]
    genome.seek(0)
    for line in genome:
        if line.startswith(">"):
            #search for chromo of key
            if dictChromo == line[1:-1]:
                read = True
            else:
                read = False
        else:
            if read == True:
                #search for coord key in genomic chromosome seq
                tempLen += len(line[:-1])
                if tempLen >= int(coordDict[key][0][0]):
                    start = int(coordDict[key][0][0]) - 1
                    stop = int(coordDict[key][0][1]) - 1
                    #this way only works for 1 line per seq fasta files!!!!
                    #seq = line[start:stop]
                    #finalDict[key] = [seq,]
                lastTempLen = tempLen
                    
                    





pseu.close()
genome.close()
out.close()
