#This script is intended to add species names to my query and database gene
#names in my .blast file.
#Created by David E. Hufnagel 01-04-2012

import sys

blast = open(sys.argv[1])    #the name fo the .blast file to be processed
out = open(sys.argv[2], "w") #the name of the output .blast file
name1 = sys.argv[3]          #the word to append to the beginning of all words in the first column
name2 = sys.argv[4]          #the word to append to the beginning of all words in the second column

#print the unix command line that called this script
out.write('#python %s\n'%(' '.join(sys.argv)))

for line in blast:
    #skip the first info line if it's there
    if not line.startswith("#"):
        lineLst = line.split("\t")
        newName1 = name1 + "_" + lineLst[0]
        newName2 = name2 + "_" + lineLst[1]
        newLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (newName1, newName2, lineLst[2], lineLst[3], lineLst[4], lineLst[5], lineLst[6], lineLst[7], lineLst[8], lineLst[9], lineLst[10], lineLst[11])
        out.write(newLine)

blast.close()
out.close()
