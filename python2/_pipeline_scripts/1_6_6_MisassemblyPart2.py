#This script is part of a series designed to test the rates of putative genomic
#missasembly. More specifically, it takes the distance between the EST match and
#the end of the EST (Y), the distance between the matching region on the genome
#and the end of the contig (X) and the average intron length for the species of
#interest (Z) to calculate whether Y <= X + Z  If true, add a missasemby
#evidence to the contig, with a bit more complexity.  Part II is designed to take
#the part1 output file with 2 lines per est and condense it into one line based
#on the following logic: 1) If either EST is labeled as "misassembled", the
#group is "misassembled" 2) If either EST is "good", but neither is "misassembled"
#the group is "good" and 3) If both ESTs are "ambiguous" the group is "ambiguous"
#the condensed info is output info in the following format:
#estName   contigName   misassembled?
#with one line per est
#Created by David E. Hufnagel on August 21, 2012

import sys

inp = open(sys.argv[1])      #input double form part 1 est match file
out = open(sys.argv[2], "w") #output condesed part 2 est match file




def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)



        
#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Import input file into dict of key: (estName, contigName) val: [isMisassembled,]
inpDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict((lineLst[0], lineLst[1]), lineLst[2].strip(), inpDict)

#Go through dict, make decisions and output info
for tup in inpDict:

    #Quality check to make sure there's only 2 lines per EST
    if len(inpDict[tup]) != 2:
        print "\n***more than 2 lines per EST***\n"
        print tup
        print inpDict[tup]

    #Go through isMisassemblys, and set new isMisassembled
    if "misassembled" in inpDict[tup]:
        isMisassembled = "misassembled"
    elif "good" in inpDict[tup]:
        isMisassembled = "good"
    elif inpDict[tup] == ["ambiguous","ambiguous"]:
        isMisassembled = "ambiguous"
    else:
        print inpDict[tup]
        print "ERROR HERE"

    estName = tup[0]
    contigName = tup[1]
    newLine = "%s\t%s\t%s\n" % (estName, contigName, isMisassembled)
    out.write(newLine)




inp.close()
out.close()
