#This script is designed to change the names in the first 2 columns in a
#blast file so that they match fasta file names
#WARNING VERY SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])       #The input m8 blast file
out = open(sys.argv[2], "w")  #The output m8 blast file


def FixName(oldName):
    newName = ""
    temp = []
    
    if oldName.startswith("Al"):
        temp = oldName.split("|")
        newName = "Aly|%s$%s$%s" % (temp[0], temp[1], temp[2])
        
    elif oldName.startswith("AT"):
        newName = "Ath|" + oldName
    
    elif oldName.startswith("Bra"):
        newName = "Bra|" + oldName
    
    elif oldName.startswith("Rr"):
        temp = oldName.split("|")
        newName = "Rra|%s$%s" % (temp[0], temp[1])
    
    else:
        print "***ERROR HERE***"
    return newName


for line in inp:
    lineLst = line.split("\t")
    
    name1 = lineLst[0]
    newName1 = FixName(name1)

    name2 = lineLst[1]
    newName2 = FixName(name2)

    theRest = "\t".join(lineLst[2:])
    newLine = "%s\t%s\t%s" % (newName1, newName2, theRest)
    out.write(newLine)



inp.close()
out.close()
