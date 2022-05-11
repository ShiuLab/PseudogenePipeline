#Designed to take .tblastn_parsed_G1451.PE_I1451.PS1, .tblastn_parsed and
#.4col.true.RMfilt.hiConf.cdnm.noTransp (final) pseudogene files and make a
#pseudoexon file in the format:
#name  pseudogene  chromo  start  stop  %ID e-val score
#Created by David E. Hufnagel on June 24, 2013
import sys

blast = open(sys.argv[1])      #input pre-pseudogene pipeline m8 parsed .tblastn_parsed file
mergedPseu = open(sys.argv[2]) #input .tblastn_parsed_G1451.PE_I1451.PS1 just merged pseudogenes file
finalPseu = open(sys.argv[3])  #input fully filtered, 4col pseudogenes file
out = open(sys.argv[4], "w")   #output pseudoexons file



def OrderCoords(coor1, coor2):
    coor1 = int(coor1);coor2 = int(coor2)
    
    if coor1 < coor2:
        small = coor1
        big = coor2
    elif coor2 < coor1:
        small = coor2
        big = coor1
    else:
        print "***ORDERCOORDS ERROR HERE***"
        return coor1,coor2

    return small, big

def GetMinMax(listx):
    #build the list
    tempLst = [[]]
    cnt = 0
    for ch in str(listx):
        if ch.isdigit():
            tempLst[cnt].append(ch)
        elif ch == ",":
            tempLst[cnt] = "".join(tempLst[cnt])
            tempLst.append([])
            cnt += 1
    else:
        tempLst[cnt] = "".join(tempLst[cnt])

    #make list items ints and return the min and max
    cnt = 0
    for item in tempLst:
        tempLst[cnt] = int(item)
        cnt += 1

    return min(tempLst), max(tempLst)

def GetPairs(listx):
    #build the list
    tempLst = [[]]
    cnt = 0
    for ch in str(listx):
        if ch.isdigit():
            tempLst[cnt].append(ch)
        elif ch == ",":
            tempLst[cnt] = "".join(tempLst[cnt])
            tempLst.append([])
            cnt += 1
    else:
        tempLst[cnt] = "".join(tempLst[cnt])

    #make list items ints and break them up into groups of 2
    cnt = 0
    for item in tempLst:
        tempLst[cnt] = int(item)
        cnt += 1

    pairs = zip(tempLst[::2], tempLst[1::2])
    return pairs

def SaveIntoDict(key, val, dictX):
    if key not in dictX:
        dictX[key] = [val]
    else:
        dictX[key].append(val)


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through finalPseu and make a dict of key: chromo;start;stop val: pgName
pgDict1 = {}
for line in finalPseu:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        pgName = lineLst[0]
        chromo = lineLst[1]
        pgStart,pgStop = OrderCoords(lineLst[2], lineLst[3])
        SaveIntoDict(chromo, (pgStart, pgStop, pgName) ,pgDict1)


#Go through mergedPseu and add to the dict so that it's key: pgName
#val: (chromo, start, stop, [(peStart1, peStop1),(peStart2, peStop2),(peStart3, peStop3)...])
#and make sure all pseudogenes are found in this file
pgDict2 = {}
for line in mergedPseu:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        #print lineLst
        #print lineLst[2]
        pgStart, pgStop = GetMinMax(lineLst[2])
        #if line.startswith("Mescscaffold11279\tcassava4_1_008010"):
            #print lineLst
            #print pgStart
            #print pgStop
            #sys.exit()
        chromo = lineLst[0]
        #print group
        if chromo in pgDict1:# and line.startswith("Mescscaffold11279\tcassava4_1_008010"):
            #print lineLst
            #print pgDict1[chromo]
            for pseudo in pgDict1[chromo]:
                #print pseudo
                #print pgStart
                #print pgStop
                #print pseudo[0]
                #print pseudo[1]
                if (pgStart >= int(pseudo[0]) and pgStart <= int(pseudo[1])) or (pgStop >= int(pseudo[0]) and pgStop <= int(pseudo[1])):
                    #print "in"
                    pgName = pseudo[2]
                    pgDict2[pgName] = (chromo, pgStart, pgStop, [])
                    #fill the list with pseudoexons
                    cnt = 0
                    for pair in GetPairs(lineLst[2]):
                        small, big = OrderCoords(pair[0], pair[1])
                        pgDict2[pgName][3].append((small, big))
                        cnt += 1
                    #print pgDict2[pgName]
            #print
        
print len(pgDict1)
print len(pgDict2)
print
for name in pgDict2:
    print name, pgDict2[name]
        

#Go through blast, find pseudoexons, gather perID, e-val and score for each
#pseudoexon and output info in the format:name  pseudogene  chromo  start  stop  %ID e-val score






blast.close()
mergedPseu.close()
finalPseu.close()
out.close()
