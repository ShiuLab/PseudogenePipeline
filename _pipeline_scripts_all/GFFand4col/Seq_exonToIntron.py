#This script is designed to make a 4 col intron file from a 4 col exon file in
#the standard 4 col format (name, chromo, start, stop)
#Created by David E. Hufnagel on June 25, 2012
"""Algorithm:
1) go through inp and collect exons into groups [a dict of lineLsts (the name
   is the key and the value is a list of lineLsts w/o the names)]
2) Determine new coordinates and Put the new coordinates into a new dict of lineLsts (the name is the key
   and the value is a list of lineLsts w/o the names)
3) Output the new lines"""


import sys

inp = open(sys.argv[1])      #the input exon file
out = open(sys.argv[2], "w") #the output intron file





def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)
        
def OrderCoords(oldCoor1, oldCoor2):
    coor1 = "";coor2 = ""
    coorA = int(oldCoor1);coorB = int(oldCoor2)

    if coorA < coorB:
        coor1 = coorA
        coor2 = coorB
    elif coorB < coorA:
        #print "###########################################################"
        coor1 = coorB
        coor2 = coorA
    else:
        print "\n***ERROR HERE***"
        coor1 = coorA
        coor2 = coorB
        print coorA
        print coorB
        print

    return coor1, coor2

def AddZeroes(num): # 3  -> "0003"
    tempNum = str(num)
    while len(tempNum)<4:
        tempNum = "0" + tempNum    
    
    newNum = tempNum
    return newNum

def SortGroup(group):
    nGroup = []
    for x in group:
        nCoor1, nCoor2 = OrderCoords(x[1], x[2].strip())
        nX = [x[0], nCoor1, nCoor2]
        nGroup.append(nX)
    nGroup.sort(key=lambda member: member[1])

    return nGroup
    


        
#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#1) go through inp and collect exons into groups [a dict of lineLsts (name is the key and the value is a list of lineLsts w/o the names)]
oldDict = {}
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")          ###NOTE THESE NEXT 2 LINES###
        name = lineLst[0].split("Exon")[0]   #general purpose
        #name = lineLst[0].split(".")[0]     #specifically for Br
        SaveIntoDict(name, lineLst[1:], oldDict)

#2) Determine new coordinates and Put the new coordinates into a new dict of
#lineLsts (the name is the key and the value is a list of lineLsts w/o the names)
for name in oldDict:
    group = oldDict[name]
    tmpCoordLst = []  #will hold a list of intron coordinates for this protein
    lastBig = 0       #the bigger value in the coordinates of the last exon
    newGroup = SortGroup(group)
    for lineLst in newGroup:
        small, big = OrderCoords(lineLst[1], lineLst[2]) #the stripping of lineLst[2] was handled in SortGroup()
        
        if lastBig != 0:
            tmpCoordLst.append((lineLst[0], lastBig + 1, small - 1))

        lastBig = big

    num = 1
    for newLst in tmpCoordLst:
        newNum = AddZeroes(num)
        newName = "%s%s%s" % (name, "Intron", newNum)
        newLine = "%s\t%s\t%s\t%s\n" % (newName, newLst[0], newLst[1], newLst[2])
        out.write(newLine)
        num += 1




inp.close()
out.close()
