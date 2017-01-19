#This script is designed to take a final 4 column pseudogene file and make a GFF
#file from it.
#Created by David E. Hufnagel on June 29, 2012

import sys

four = open(sys.argv[1])     #input 4col pseudogene file with code-derived name      (.real4col.RMfilt.cdnm.noCM)
ref = open(sys.argv[2])      #reference file with original long-form pseudogene name (.real4col.mod.fa.ref)
out = open(sys.argv[3], "w") #output GFF file                                        (.GFF)





def ForwardOrReverse(start, stop):
    if int(start) < int(stop):
        return "+"
    elif int(start) > int(stop):
        return "-"
    else:
        return "***FR error!!!***"

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)




        
#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through ref and put it into a dict
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[0], lineLst[1].strip(), refDict)

#Go through four and output the file into a GFF format
for line in four:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        name = lineLst[0]
        chromo = lineLst[1]
        source = "ShiuLab"
        typex = "pseudogene"
        start = lineLst[2]
        stop = lineLst[3].strip()
        score = "."
        strand = ForwardOrReverse(start, stop)
        phase = "."
        code = "".join(refDict[name[2:].split("_")[0]][0].split(";")[-1].split("|"))
        attributes = "ID=%s;Name=%s;Note=evidence code is %s. "\
                     "functional protein paralog is the second part of the "\
                     "pseudogene name." % (name, name, code)
        newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chromo, source, typex, start, stop, score, strand, phase, attributes)
        out.write(newLine)




four.close()
ref.close()
out.close()
