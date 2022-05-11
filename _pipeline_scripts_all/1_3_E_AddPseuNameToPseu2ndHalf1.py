#This script is designed to take the goofy names in the .ps_cds.2nd_half files
#and replace them with pseudogene names. Input cds names are in the
#format: >Alyrscaffold_8|6837211-6837357_1|147
#Created by David E. Hufnagel on Mar 21, 2013
import sys

cds = open(sys.argv[1])
ref = open(sys.argv[2])
out = open(sys.argv[3], "w")



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and make a dict of key: goofyName val: codeName
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        codeName = lineLst[0]
        bigName = lineLst[1]
        goofyName = bigName.split(";")[1]+ "_" + bigName.split(";")[2].split(":")[1].replace("-","|")
        refDict[goofyName] = codeName
print refDict

#Go through cds, do the translation and output result
for line in cds:
    if not line.startswith("#"):
        if line.startswith(">"):
            currName = line.strip(">").strip("\n")
            print codeName
            codeName = refDict[currName]
            line = ">%s\n" % (codeName)
        out.write(line)




cds.close()
ref.close()
out.close()
