#This script is designed to compare Brapa_197_annotation_info.txt to
#at_br_block.csv
#Created by David E. Hufnagel on 4-10-2012
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

phyt = open(sys.argv[1])     #Brapa_197_annotation_info.txt
csv = open(sys.argv[2])      #at_br_block.csv
out = open(sys.argv[3], "w") #output file


#go through phyt
phytDict = {}
for line in phyt:
    lineLst = line.split("\t")
    if lineLst[6] != "":
        if lineLst[0] not in phytDict:
            phytDict[lineLst[0]] = lineLst[6][:-2]
        else:
            print "Heeyyyyy!!!!!"  #all 1:1!

#go through csv
csvDict = {}
csv.readline()
for ln in csv:
    lnLst = ln.split(",")
    if lnLst[4] not in csvDict:
        csvDict[lnLst[4]] = [lnLst[3]]
    else:
        csvDict[lnLst[4]].append(lnLst[3])

###compare phyt to csv
##yes = 0
##no = 0
##for l in phytDict:
##    if l in csvDict:
##        if phytDict[l] in csvDict[l]:
##            yes += 1
##        else:
##            no += 1
##    else:
##        no += 1

###compare csv to phyt
##yes = 0
##no = 0
##for l in csvDict:
##    for x in csvDict[l]:
##        if l in phytDict:
##            if x == phytDict[l]:
##                yes += 1
##            else:
##                no += 1
##        else:
##            no += 1
##    else:
##        print

print yes
print no
    
    









    



phyt.close()
csv.close()
out.close()
