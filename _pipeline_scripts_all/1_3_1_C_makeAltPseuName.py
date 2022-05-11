#This script is designed to make a pseudogene fasta file with an alternative
#name that works with Shin-Han's script_step7.py and script_step7.1.py
#Created by David E. Hufnagel on Jan 14, 2013
#Updated on Mar 19, 2013 to handle proper .ref files

import sys

ref = open(sys.argv[1])      #input reference file connecting original names to big names (alternative names will be developed from the big names)
fasta = open(sys.argv[2])    #input fasta file with code names
out = open(sys.argv[3], "w") #output fasta file with alternative names



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#import ref into a dict of key: original name val: alt name
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        FG = lineLst[1].split(";")[0]#.split("|")[0] #the extra stuff is for Rr only
        #fullCodeName = "Ps%s_%s" % (lineLst[0], FG)
        fullCodeName = lineLst[0]
        altName = lineLst[1].split(";")[1] + "|" + lineLst[1].split(";")[2].split(":")[1]
        refDict[fullCodeName] = altName

#go through fasta and output lines with alternative names
for line in fasta:
    if not line.startswith("#"):
        if line.startswith(">"):
            origName = line.strip(">").strip("\n")
            altName = refDict[origName]
            out.write(">%s\n" % (altName))
        else:
            out.write(line)






ref.close()
fasta.close()
out.close()
