#This script is designed to take a 4 col file with a real name in the first
#column and replace it with a derivitive of code name that's already been
#generated.  Specifically for making new pseudogene names for Gaurav to handle.
#Created by David E. Hufnagel on June 27, 2012
#WARNING: SOMEWHAT SITUATION SPECIFIC

import sys

four = open(sys.argv[1])     #the 4 col file with the pseudogene info and the real names
ref = open(sys.argv[2])      #the reference file with the real names and code names
out = open(sys.argv[3], "w") #ouput 4 col file with code-derived names


def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = val
    else:
        print "error"


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and put it into a dict
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[1].strip(), lineLst[0], refDict)
        
#Go through four and make new file
for line in four:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        prot = lineLst[0].split(";")[0].split("|")[0]
        name = "Ps%s_%s" % (refDict[lineLst[0]], prot)
        newLine = "%s\t%s\t%s\t%s" % (name, lineLst[1], lineLst[2], lineLst[3])
        out.write(newLine)






four.close()
ref.close()
out.close()
