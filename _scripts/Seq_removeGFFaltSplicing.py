#This script is designed to take a gff file and remove all alternative splice
#mRNAs, UTRs and CDSs
#Created by David E. Hufnagel on Feb 27, 2013
import sys

gff = open(sys.argv[1])      #input gff file
out = open(sys.argv[2], "w") #output gff file without alternative splicing



#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through gff and make a list of PACIDs for longest transcripts
goodLst = []
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "mRNA":
            isLongest = bool(int(lineLst[8].split("longest=")[1].split(";")[0]))
            pacID = lineLst[8].split("pacid=")[1].split(";")[0]
            if isLongest == True:
                goodLst.append(pacID)

#Go through gff again to output longest transcripts
gff.seek(0)
for line in gff:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "gene":
            out.write(line)
        elif lineLst[2] == "exon" or lineLst[2] == "mRNA" or lineLst[2] == "CDS" or "UTR" in lineLst[2]:
            pacID = lineLst[8].split("pacid=")[1].split(";")[0]
            if pacID in goodLst:
                out.write(line)



gff.close()
out.close()
