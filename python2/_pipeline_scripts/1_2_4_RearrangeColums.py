#This script is designed to shift columns around in a tab delimited file
#Created by: David E. Hufnagel (10-13-2011)
#Updated on 4-17-2012: allowed for conversions if different sizes (Ex: instead
#of only 1,2,3 -> 2,1,3 can do 1,2,3,4 -> 3,1,2)

import sys

OLD = sys.argv[1]  #old filename
orderStr = sys.argv[2] #new order of columns.  Comma delimited.  Ex: 3,1,2,4,5
NEW = sys.argv[3]  #new filename
old = open(OLD)
new = open(NEW, "w")

#goes through a list and subtracts 1 from each integer value
def NumToIndex(listx):
    cnt = 0
    for n in listx:
        listx[cnt] = int(n)-1
        cnt += 1
    return listx

#writes the users command line prompt on the first line of the output file.
new.write("#python %s\n" % (" ".join(sys.argv)))

#turn number list into python index values (subtract 1)
orderLst = orderStr.split(",")
orderLst = NumToIndex(orderLst)

for line in old:
    if not line.startswith("#"):
        tempLst = []
        lineLst = line.split("\t")
        
        for n in orderLst:
            tempLst.append(lineLst[n].strip("\n"))

            #write line in new order
            temp = "\t".join(tempLst)
            #print temp
        new.write(temp + "\n")
    
old.close()
new.close()
