#This script is designed to make a pseudogene fasta file with the proper code
#names from the alternative name developed for Shin-Han's script_step7.py and
#script_step7.1.py.
#Created by David E. Hufnagel on Jan 14, 2013
#Updated on Mar 19, 2013 to handle proper .ref files

import sys

ref = open(sys.argv[1])      #input reference file connecting original names to big names (alternative names will be developed from the big names)
fasta = open(sys.argv[2])    #input fasta file with alternative names
out = open(sys.argv[3], "w") #output fasta file with original names

#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Import ref into a dict of key: alt name val: orig name
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        FG = lineLst[1].split(";")[0]
        #fullCodeName = "Ps%s_%s" % (lineLst[0], FG)
        fullCodeName = lineLst[0]
        altName = lineLst[1].split(";")[1] + "|" + lineLst[1].split(";")[2].split(":")[1]
        refDict[altName] = fullCodeName

#Go through fasta file and output the file with original (rather than alternative) names
for line in fasta:
    if not line.startswith("#"):
        if line.startswith(">"):
            #print altName
            altName = line.strip(">").strip("\n")
            altName = "%s|%s" % ("_".join(altName.split("_")[:-1]), altName.split("_")[-1])
            altName = altName.split("|")
            altName = "|".join(altName[:-1]) + "-" + altName[-1]
            #print altName
            origName = refDict[altName]#.split("|")[0] #the extra stuff is for Rr only
            #print origName
            out.write(">%s\n" % (origName))
        else:
            out.write(line)




ref.close()
fasta.close()
out.close()
