#This script is designed to make a 4col proteins file from a masked genome file
#where protein coding genes are masked as "M"s
#Created by David E. Hufnagel on Aug 7, 2012

import sys

fasta = open(sys.argv[1])     #input masked fasta file
four = open(sys.argv[2], "w") #output 4col proteins file derived from masked regions





def MainWork(cnt, chromo):
    mStart = seq.find("M")
    mStop = seq.rfind("M")
    mSeq = seq[mStart:mStop+1]
                
    #no proteins
    if mSeq == "":
        pass #this is necessary because when Mseq = "" mSeq = len(mSeq) * "M"
        
    #one protein
    elif mSeq == len(mSeq) * "M":
        fakeName = "AtFakeProt%s" % (str(cnt).zfill(6))
        newLine = "%s\t%s\t%s\t%s\n" % (fakeName, chromo, mStart + 1, mStop + 1)
        four.write(newLine)
        cnt += 1
        
    #multiple proteins
    else:
        coordCnt = mStart #the count to keep track of the position on the contig while going through the seq
        mCoordLst = []
        tempStart = mStart
        inMs = True
        for base in mSeq:
            #if coordCnt == mStart: #on the first iteration, set the first tuple
            if base != "M":
                if inMs == True:
                    #find the first stop and add forst tuple
                    mCoordLst.append((tempStart, coordCnt - 1))
                inMs = False
            else:
                if inMs == False:
                    tempStart = coordCnt
                inMs = True
            coordCnt += 1
        #on the way out add the last tuple and write into the output
        else:
            mCoordLst.append((tempStart, coordCnt - 1))
            for tup in mCoordLst:
                fakeName = "AtFakeProt%s" % (str(cnt).zfill(6))
                newLine = "%s\t%s\t%s\t%s\n" % (fakeName, chromo, tup[0] + 1, tup[1] + 1)
                four.write(newLine)
                cnt += 1

    return cnt





chromo = ""   #the current name associated with a seq
seq = ""        #the current seq to be built up for each seq line
cnt = 1         #the count for use in making fake protein names
for line in fasta:
    if not line.startswith("#"):
        if line.startswith(">"):
            if chromo != "":                
                #see if the whole thing is Ms and handle them (AKA there's only one protein coding region there)
                #no proteins
                cnt = MainWork(cnt, chromo)

            chromo = line.strip().strip(">")
            seq = ""
            
        #add up seq on each line until it's the whole thing
        else:
            seq += line.strip()
            
#handle the last seq on the way out
else:
    cnt = MainWork(cnt, chromo)



fasta.close()
four.close()
