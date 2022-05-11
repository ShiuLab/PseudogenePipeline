#This script is designed to take a gff3 file and from it make a gff of CDSs    
#that are part of the longest transcripts in the gene
#Created by David E. Hufnagel on May 2, 2013
import sys

inp = open(sys.argv[1])      #input gff3 file
out = open(sys.argv[2], "w") #output gff file of primary transcript CDSs


#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Do the actual processing
longDict = {} #a dict of key: PACID val: isLongest for each mRNA
for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "gene":
            pass
        #make a dict of key: PACID val: isLongest
        elif lineLst[2] == "mRNA":
            pacid = lineLst[8].split("pacid=")[1].split(";")[0]
            isLongest = bool(int(lineLst[8].split("longest=")[1].split(";")[0]))
            if pacid in longDict:
                print "ERROR: duplicate PACID"
            else:
                longDict[pacid] = isLongest
        #check if CDS is part of the primary transcript
        elif lineLst[2] == "CDS":
            #if it is output the line
            pacid = lineLst[8].split("pacid=")[1].split(";")[0]
            isLongest = longDict[pacid] 
            if isLongest:
                out.write(line)



inp.close()
out.close()
