#This script is designed to add a 5th column to my pseudogene real4col files,
#which is the second column (BLAST match subject info) in the disable_count file
#Created by David E. Hufnagel on July 11, 2012
"""Algorithm:
1) Go through reference file and save it into refDict
2) Go through dis and
    a) get code name from big name using refDict
    b) save needed info into a dict (key: codeName val: new 5th col line)
3) Go through four and add the last line using dict from part 2 (AKA output result)"""

import sys

four = open(sys.argv[1])     #input 4 col pseudogene file                   (.real4col.RMfilt.cdnm.noCM)
dis = open(sys.argv[2])      #input disable_count pseudogene file           (.fullyFiltered.disable_count.RMfilt)
ref = open(sys.argv[3])      #reference file with long names and code names (.real4col.mod.fa.ref)
out = open(sys.argv[4], "w") #output 5 col pseudogene file                  (.5col)




def FixProt(oldName):
    new = oldName.split(";")[0].split("|")[0]
    newName = "%s;%s" % (new, ";".join(oldName.split(";")[1:]))
    return newName




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) Go through reference file and save it into refDict
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        codeName = "Ps%s_%s" % (lineLst[0], lineLst[1].split(";")[0])
        bigName = lineLst[1].strip()
        bigName = FixProt(bigName)
        refDict[bigName] = codeName.split("|")[0]
        
#2) Go through dis and
#    a) get code name from big name using refDict
#    b) save needed info into a dict (key: codeName val: new 5th col line)
newColDict = {}
dis.readline()
for line in dis:
    if not line.startswith("#p"):
        if line.startswith("#"):
            lineLst = line.split(" ")
            bigName = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:])[:-1])[1:]
            bigName = FixProt(bigName)
            codeName = refDict[bigName]
            new5col = lineLst[1]
            newColDict[codeName] = new5col

#3) Go through four and add the last line using dict from part 2 (AKA output result)
for line in four:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        newLine = "%s\t%s\t%s\t%s\t%s\n" % (lineLst[0], lineLst[1], lineLst[2], lineLst[3].strip(), newColDict[lineLst[0]])
        out.write(newLine)




four.close()
dis.close()
ref.close()
out.close()
