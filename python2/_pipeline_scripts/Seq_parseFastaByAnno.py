#This script is designed to parse a cDNA fasta file by a keyword in a 4 col
#file with annotation information.
#Created by David E. Hufnagel on June 8, 2012

import sys

anno = open(sys.argv[1])     #4col file with annotation info
inp = open(sys.argv[2])      #input unfiltered fasta file
out = open(sys.argv[3], "w") #output filtered fasta file
key = sys.argv[4]            #key word for what to keep in the filtration process





#go through anno and make dict of key:name value:annotation
annoDict = {}
for line in anno:
    lineLst = line.split("\t")
    name = lineLst[0].split("|")[0]
    print name
    annot = lineLst[0].split("|")[1]
    annoDict[name] = annot

#go through inp and filter output by key
write = True
for line in inp:
    if line.startswith(">"):
        try:
            if annoDict[line[1:].strip()] == key:
                write = True
                out.write(line)
            else:
                write = False
        except KeyError:
            print "\nhandled KeyError\n: ", line[1:].strip()
    else:
        if write == True:
            out.write(line)






anno.close()
inp.close()
out.close()
