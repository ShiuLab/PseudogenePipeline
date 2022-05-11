#This script is designed to get the number of exons/introns per gene.  The
#algorithm is simle, counting the number of different items in column 2 per
#unique item in column 1.
#Created by David E. Hufnagel on July 3, 2012

import sys

inp = open(sys.argv[1])      #input exon/intron file
out = open(sys.argv[2], "w") #output data file
ID = sys.argv[3]             #two char species ID followed by In or Ex:  AtIn





def SaveIntoNumDict(key, dictX):
    if key not in dictX:
        dictX[key] = 1
    else:
        dictX[key] += 1





if len(ID) != 4:
    print "ID is incorrect!"
    
#Writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Save input exon/intron file proteins into a number dictionary 
inpDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if ID.endswith("Ex"):
            if ID.startswith("Br"):
                SaveIntoNumDict(lineLst[0].split(".CDS.")[0], inpDict)
            else:
                SaveIntoNumDict(lineLst[0].split("Exon")[0], inpDict)
        elif ID.endswith("In"):
            SaveIntoNumDict(lineLst[0].split("Intron")[0], inpDict)
        else:
            print "\n***error here***\n"

#Output data into a file
for key in inpDict:
    val = inpDict[key]
    newLine = "%s\t%s\n" % (key, val)
    out.write(newLine)





inp.close()
out.close()
