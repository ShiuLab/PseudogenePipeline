#This script was designed to organize and improve Gaurav's gff file for R.
#raphanistrum.  This task includes 1) filtering out useless crap (according to
#columns two and three) 2) changing the names to match all our other data 3)
#organize the file so that it looks like a normal gff instead of all genes being
#clustered together and all exons clustered together and etc.  4) Calculate the
#longest isoform and label it with the tag longest=1 and label the others with
#the tag longest=0
#Created by David E. Hufnagel on March 6, 2013
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

gff = open(sys.argv[1])
ref = open(sys.argv[2])
out = open(sys.argv[3], "w")



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and make a dict of key: oldName val: newName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        refDict[lineLst[0]] = lineLst[2]

#Go through gff
tempSet = set()
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        tempSet.add("%s___%s" % (lineLst[1], lineLst[2]))
        #?
for item in tempSet:
    print item




gff.close()
ref.close()
out.close()
