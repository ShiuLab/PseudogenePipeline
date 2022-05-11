#This script is designed to make a 4 column file of the ends of contigs given
#the fasta file for the contigs
#Created by David E. Hufnagel on June 12, 2012
"""Algorithm:
1) Go through the input file
    a) makeup name
    b) gather chromo
    c) get coords
    d) output 4col line"""

import sys

fasta = open(sys.argv[1])    #the input fasta file (contigs)
out = open(sys.argv[2], "w") #the output 4col file with the ends of the contigs
dist = int(sys.argv[3])      #the distance from the end to be concerned with





def ImportFastaToDict(fasta):
    name = ""
    seq = ""
    fastDict = {}
    for line in fasta:
        if line.startswith(">"):
            fastDict[name] = seq
            name = line[1:-1]
            seq = ""
        else:
            seq += line[:-1]
    #for the last name, seq pair
    else:
        fastDict[name] = seq
        fastDict.pop("")

    return fastDict
            





#Go though fasta file and inport it into a dict with key: name value: seq
fastDict = ImportFastaToDict(fasta)

start1 = "";stop1 = ""
start2 = "";stop2 = ""
for key in fastDict:
    val = fastDict[key]
    if len(val) >= dist:
        start1 = 1
        stop1 = dist
        name1 = "%s;%s-%s" % (key, start1, stop1)
        
        start2 = len(val) - dist + 1
        stop2 = len(val)
        name2 = "%s;%s-%s" % (key, start2, stop2)
        
        begLine = "%s\t%s\t%s\t%s\n" % (name1, key, start1, stop1)
        endLine = "%s\t%s\t%s\t%s\n" % (name2, key, start2, stop2)
        out.write(begLine)
        out.write(endLine)
    else:
        start2 = start1 = 1
        stop2 = stop1 = len(val)
        name2 = name1 = "%s;%s-%s" % (key, start1, stop1)
        
        begLine = "%s\t%s\t%s\t%s\n" % (name1, key, start1, stop1)
        endLine = "%s\t%s\t%s\t%s\n" % (name2, key, start2, stop2)
        out.write(begLine)
        out.write(endLine)







fasta.close()
out.close()

