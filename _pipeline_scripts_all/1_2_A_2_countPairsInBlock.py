#This script is designed to count the number of gene pairs in each syntenic
#block found in an MCScanX .collinearity file and export the lengths into
#a tab delimited file that could be imported in excel for generating a
#histogram.  The algorithm used is that it reads the information from the input
#file and puts information for the number of genes in a block in a dictionary
#where the key is the number of genes and the value is the number of blocks
#fitting that description.
#Created by David E. Hufnagel on 3-12-2012

import sys

inp = open(sys.argv[1])      #input .collinearity file
out = open(sys.argv[2], "w") #output file



#skip past the title in the .collinearity file
for ln in inp:
    if not ln.startswith("#"):
        break

#extract the info from the .collinearity file
bigDict = {}
pairCnt = 0
name = ""
lastName = ""
for line in inp:
    if line.startswith("##"):  #block info line
        name = line.split(":")[0][13:]
        if pairCnt != 0:       #skipping the first block because data is gathered one block behind
            if int(pairCnt) not in bigDict:
                bigDict[int(pairCnt)] = 1
            else:
                bigDict[int(pairCnt)] += 1
            pairCnt = 0   
        lastName = name
    else:                      #gene pair info line
        pairCnt += 1
        
else:                          #end of file
    if int(pairCnt) not in bigDict:
        bigDict[int(pairCnt)] = 1
    else:
        bigDict[int(pairCnt)] += 1        

#export the info in a tab delimited text file
for key in bigDict:
    newLine = "%s\t%s\n" % (key, bigDict[key])
    out.write(newLine)





inp.close()
out.close()
