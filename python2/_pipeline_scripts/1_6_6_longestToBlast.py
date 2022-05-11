#This script is designed to convert a 4col .longest parsed GMAP file to a
#tabular m8 blast format.  This requires the .breaks file
#Created by David E. Hufnagel on Sep 12, 2012

import sys

breaks = open(sys.argv[1])   #input .breaks partially parsed GMAP output file
longest = open(sys.argv[2])  #input .longest 4col parsed GMAP output file
out = open(sys.argv[3], "w") #output blast format parsed GMAP output file




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through breaks file and put file into a list of lists
breakLst = []
tempLst = []
for line in breaks:
    if not line.startswith("#Name"):
        if line.startswith("##"):
            if not tempLst == []:
                breakLst.append(tempLst)
                tempLst = []
        else:
            tempLst.append(line)

#go through breaks list and gather query match information into a dict of
#key: query;database value: startQueryCoord-endQueryCoord
breakDict = {}
for group in breakLst:
    #if there's only one member of the group just take the value of the match
    lineLst = group[0].split("\t")
    key = lineLst[0] + ":" + lineLst[7].split(":")[0]
    if len(group) == 1:
        queryCoords = lineLst[5]
        breakDict[key] = queryCoords

    #if there are many members generate a minmax coord
    else:
        queryStart = 9999999999
        queryStop = 0
        for line in group:
            lineLst = line.split("\t")
            tempStart = int(lineLst[5].split('-')[0])
            tempStop = int(lineLst[5].split('-')[1])
            if tempStart < queryStart:
                queryStart = tempStart
            if tempStop > queryStop:
                queryStop = tempStop
            queryCoords = "%s-%s" % (queryStart, queryStop)
            breakDict[key] = queryCoords

#go through longest file and make the output blast format parsed GMAP output file from the input information
for line in longest:
    if not line.startswith("#"):
        lineLst = line.split()
        query = lineLst[0]
        subj = lineLst[1]
        key = query + ":" + subj
        qStart = breakDict[key].split("-")[0]
        qEnd = breakDict[key].split("-")[1]
        length = abs(int(qEnd) - int(qStart))
        sStart = lineLst[2]
        sEnd = lineLst[3]
        newLine = "%s\t%s\t-\t%s\t-\t-\t%s\t%s\t%s\t%s\t-\t-\n" % (query, subj, length, qStart, qEnd, sStart, sEnd)
        out.write(newLine)



breaks.close()
longest.close()
out.close()
