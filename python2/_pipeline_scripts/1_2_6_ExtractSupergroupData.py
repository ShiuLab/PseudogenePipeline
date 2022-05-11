#This script is designed to take chunks of a CIRCOS links file seperated
#by #####whatever\n and make a CSV file for use in excell to make histograms
#on the distances between links in a superblock.
#Made by: David E. Hufnagel on 1-3-2012

import sys
import csv

inp = open(sys.argv[1])       #The input derived from a links file
out = open(sys.argv[2], "w")  #The CSV output file
name1 = sys.argv[3]           #The species name for the odd rows   Ex: Ath
name2 = sys.argv[4]           #The species name for the even rows  Ex: Aly
writer = csv.writer(out, dialect = "excel")

def ProcessRows(line, cnt):
    lineLst = line.split("\t")

    #odd rows (first species)
    if cnt%2 != 0:
        #print lineLst[1]    ####always uncomment this at first to make sure it is getting the right species throughout
        list1.append(lineLst)

    #even rows (second species)
    elif cnt%2 == 0:
        list2.append(lineLst)
        
    else:
        print "ERROR 0!"
        
    cnt += 1
    return cnt

    
    
    

cnt = 1
currName = ""
#the list of lists that is the organized rows Structure Lv1(biggest): [supergroup1, supergroup2, supergroup3]  Lv2: [Species1, Species2]  Lv3: [row1, row2, row3...]  Lv4: [word1, word2, word3...]
bigLst = []
list1 = []  #a list of lists for species 1
list2 = []  #a list of lists for species 2
names = []  #The names of the supergroups, taken from whatever follows the "#####" on the same line.
            # example name: "Aly3_Ath6", meaning for A. Lyrata it will be supergroup 3 and for A. Thaliana it will be group 6.
for line in inp:
    if line.startswith{"#") and not line.startswith
    
    if line.startswith("#####"):
        currName = line[5:-1]
        names.append(currName)
        #reset count for determining odd and even
        cnt = 1
        
        bigLst.append([list1, list2])
        list1 = []
        list2 = []
    else:
        cnt = ProcessRows(line, cnt)

numMembers = 45
row1 = ["Summary Info:","Supergroup Name", name1, "Number of Members", numMembers]
row2 = ["Start", "Stop", "Distance"]
writer.writerow(row1)
writer.writerow(row2)

print bigLst[0]

inp.close()
out.close()
