#This script is designed for going through a karyotype file and fixing the
#chromosome names on the band lines so that they match the chromosome names
#on the chromosome lines

import sys

def fixScafName(oldName): #Al_scaffold_3  ->  Al_scaf_0003
    tempLst = oldName.split("_")

    #change number
    tempNum = tempLst[2]
    while len(tempNum)<4:
        tempNum = "0" + tempNum    
    
    newName = tempLst[0] + "_scaf_" + tempNum
    return newName

INP = sys.argv[1] #The name of the karyotype file to be filtered
OUT = sys.argv[2] #The name of the new karyotype file (this programs output)

inp = open(INP)
out = open(OUT, "w")

#print the unix command line that called this script
print out.write('#python %s\n'%(' '.join(sys.argv)))

cnt = 0
for line in inp:
    lineLst = line.split("\t")
    
    if lineLst[0] == "chr":   #chromosome lines
        out.write(line)
        
    elif lineLst[0] == "band" and lineLst[1].startswith("Al"):  #band lines
        oldName = lineLst[1]
        newName = fixScafName(oldName)
        temp = lineLst
        
        lineLst = []
        lineLst.append(temp[0])
        lineLst.append(newName)
        lineLst.append(temp[2])
        lineLst.append(temp[3])
        lineLst.append(temp[4])
        lineLst.append(temp[5])
        lineLst.append(temp[6])

        line = "\t".join(lineLst)
        out.write(line)

    else:
        out.write(line)
        
            
    cnt += 1

inp.close()
out.close()
