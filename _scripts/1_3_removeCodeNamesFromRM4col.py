#This script is designed to remove the code names from a generated 4col parsed
#repeat masker output file.
#Designed by David E. Hufnagel on June 8, 2012
#Updated on June 19, 2012 to include functionality for simple tabular formats
#like .4col, .bed, .blast etc and condiderations for hash lines

import sys

four = open(sys.argv[1])     #4input coded file                                    
ref = open(sys.argv[2])      #reference file with gene names and code names        
out = open(sys.argv[3], "w") #output 4col file with real gene names                
style = sys.argv[4]          #the file format of the input/output current options: (RM4col, simpTab)
col = int(sys.argv[5]) - 1   #the index of the column with the name to be replaced.  not used for style: RM4col

def SaveIntoDict(gene1, gene2, dictX):
    if gene1 not in dictX:
        dictX[gene1] = [gene2]
    else:
        dictX[gene1].append(gene2)

#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#save ref into a dict
refDict = {}
for line in ref:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        SaveIntoDict(lineLst[0], lineLst[1].strip(), refDict)

#go through file and replace codes with real names (output replacement)
for line in four:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if style == "RM4col":
            code = lineLst[1]
            name = refDict[code][0]
            new0 = "%s|%s" % (name, "|".join(lineLst[0].split("|")[1:]))
            newLine = "%s\t%s\t%s\t%s\t%s" % (new0, name, lineLst[2], lineLst[3], lineLst[4])
        elif style == "simpTab":
            #problem: very poorly designed method, 1) need to write a new line for each col 2) the name cannot be the last thing in a file
            code = lineLst[col]
            name = refDict[code][0]
            if col == 0:
                newLine = "%s\t%s" % (name, "\t".join(lineLst[1:]))
            elif col == 1:
                newLine = "%s\t%s" % (lineLst[0], name, "\t".join(lineLst[2:]))
            elif col == 2:
                newLine = "%s\t%s" % (lineLst[:2], name, "\t".join(lineLst[3:]))
            elif col == 3:
                newLine = "%s\t%s" % (lineLst[:3], name, "\t".join(lineLst[4:]))
            elif col == 4:
                newLine = "%s\t%s" % (lineLst[:4], name, "\t".join(lineLst[5:]))
            elif col == 5:
                newLine = "%s\t%s" % (lineLst[:5], name, "\t".join(lineLst[6:]))
            elif col == 6:
                newLine = "%s\t%s" % (lineLst[:6], name, "\t".join(lineLst[7:]))
            elif col == 7:
                newLine = "%s\t%s" % (lineLst[:7], name, "\t".join(lineLst[8:]))
            elif col == 8:
                newLine = "%s\t%s" % (lineLst[:8], name, "\t".join(lineLst[9:]))
            elif col == 9:
                newLine = "%s\t%s" % (lineLst[:9], name, "\t".join(lineLst[10:]))
            else:
                print("error:  need more col capability")
        else:
            print("***style error***")
        out.write(newLine)

four.close()
ref.close()
out.close()

print("  4-col file with gene names:", sys.argv[3])