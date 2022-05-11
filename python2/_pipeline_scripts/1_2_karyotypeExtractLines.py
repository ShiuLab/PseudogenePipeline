#This script is designed for going through a karyotype file and removing all
#but the specified bands.

import sys

INP = sys.argv[1] #The name of the karyotype file to be filtered
OUT = sys.argv[2] #The name of the new karyotype file (this programs output)
good = sys.argv[3] #the names scaffolds of the bands to be kept.  This has a
                   #specific expectation for format.  It should have no spaces
                   #and be seperated by commas. Ex: Al_scaffold_1,Al_scaffold_7

inp = open(INP)
out = open(OUT, "w")

#print the unix command line that called this script
print out.write('#python %s\n'%(' '.join(sys.argv)))

#iterate through the input file and write the acceptable lines in the
#output file.  Acceptable means a chromosome line or a band line with
#the specified name


goodLst = good.split(",")
print "your name list:"
print
print goodLst
print

cnt = 0
for line in inp:
    lineLst = line.split("\t")
    
    if lineLst[0] == "chr":   #chromosome lines
        out.write(line)
        
    if lineLst[0] == "band":  #band lines
        #write good lines
        if lineLst[1] in goodLst:
            out.write(line)
            
    cnt += 1

inp.close()
out.close()
