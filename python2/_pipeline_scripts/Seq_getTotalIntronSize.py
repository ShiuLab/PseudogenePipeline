#This script is designed to get a size file for all introns in a protein from
#the regular intron .size file
#Created by David E. Hufnagel on June 26, 2012

import sys

inp = open(sys.argv[1])      #the input 4col file
out = open(sys.argv[2], "w") #the output .size file





def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)




        
#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Make dict of lines where lines with the same name are multiple values with 1 key
bigDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        name = lineLst[0].split("Intron")[0]
        num = int(lineLst[1].strip())
        SaveIntoDict(name, num, bigDict)

#Go through the dict and calculate the total of the values. Output the new lines
for key in bigDict:
    val = bigDict[key]
    size = sum(val)
    newLine = "%s\t%d\n" % (key, size)
    out.write(newLine)





inp.close()
out.close()
