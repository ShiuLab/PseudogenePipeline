#This script is designed to add a word and a number to the string in the first
#column of a tabular file (originally designed to turn protein names into exon
#names)
#Created by David E. Hufnagel on June 19, 2012

import sys

inp = open(sys.argv[1])      #input file
out = open(sys.argv[2], "w") #output file
word = sys.argv[3]           #the word that will go in between the original name and the numbers added





def AddZeroes(num): # 3  -> "0003"
    tempNum = str(num)
    while len(tempNum)<4:
        tempNum = "0" + tempNum    
    
    newNum = tempNum
    return newNum





#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

inpDict = {} #used to count the number of times something has shown up in the file thus far
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if not lineLst[0] in inpDict:
            inpDict[lineLst[0]] = 1
        else:
            inpDict[lineLst[0]] += 1
            
        num = AddZeroes(inpDict[lineLst[0]])
        newName = "%s%s%s" % (lineLst[0], word, num)
        newLine = "%s\t%s\t%s\t%s" % (newName, lineLst[1], lineLst[2], lineLst[3])
        out.write(newLine)
    


                  
inp.close()
out.close()
