#This script is designed for converting GFF files to bed files
#Created by: David E. Hufnagel (10-13-2011)
#Updated: 12-19-2011: Added funtionality for many keywords
#(used to be just "gene")
#Updated: 3-8-2012: Changed output format to MSScanX format, made sys
#  variables more readable and added spe variable
#Updated: 6-19-2012: Made hash skipping automatic, commented out use of
#  ReverseCoords and made a new version of GetName to serve a more general
#  purpose by taking the "ID" intead of the "name".  Also switched output from
#  MCScanX standard to 4col standard

import sys

gff = open(sys.argv[1])      #GFF filename
bed = open(sys.argv[2], "w") #bed filename
keyword = sys.argv[3]        #keyword (will only keep rows with this descriptor at index 2)
spe = sys.argv[4]            #2-letter species identifier   

#print the unix command line that called this script
bed.write('#python %s\n'%(' '.join(sys.argv)))



def ReverseCoords(lineLst, coor1, coor2):
    coor1 = int(lineLst[3])
    coor2 = int(lineLst[4])
##    coorA = int(lineLst[3])
##    coorB = int(lineLst[4])
##
##    if coorA < coorB:
##        coor1 = str(coorA)
##        coor2 = str(coorB)
##    elif coorB < coorA:
##        coor1 = str(coorB)
##        coor2 = str(coorA)
##    else:
##        print "***ERROR HERE***"

    return coor1, coor2

def GetName(bunch, spe):
    ###WARNING: LOOK AT THIS BEFORE RUNNING###
    getNext = False
    if spe == "Rr":
        ID = bunch.split(";")[0].split("=")[1].split(":")[0]
    elif spe == "Br":
        ID = bunch.split(";")[1].split("=")[1]
    elif spe == "Al":
        ID = bunch.split(";")[0].split('"')[1]
    elif spe == "At":
        ID = bunch.split(";")[0].split("=")[1].strip()
    
    return ID
##    name = ""
##    chunks = bunch.split("=")
##    for chunk in chunks:
##        if chunk.endswith("Name"):
##            getNext = True
##        elif getNext == True:
##            name = chunk.strip()
##            break
##    return name
    
for line in gff:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if lineLst[2] == keyword:
            coor1 = ""
            coor2 = ""
            coor1, coor2 = ReverseCoords(lineLst, coor1, coor2)
            name = GetName(lineLst[8], spe)
            #newLine = "%s\t%s\t%s\t%s\n" % (spe + lineLst[0], name, coor1, coor2) #MCScanX
            newLine = "%s\t%s\t%s\t%s\n" % (name, spe + lineLst[0], coor1, coor2) #4col
            bed.write(newLine)

gff.close()
bed.close()
