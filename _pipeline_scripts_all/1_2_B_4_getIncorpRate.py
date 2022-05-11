#This script was designed to determine the number of Rr orthologs incorperated
#into a At-Al-Br ortholog file
#Designed by David E. Hufnagel on 5-9-2012

import sys

inp = open(sys.argv[1])  #input .comp file
orth = open(sys.argv[2]) #ortholog file, just for determining the number of lines in it



#determine number of lines in orth
orthCnt = 0
for ln in orth:
    if not ln.startswith("#"):
        orthCnt += 1

#determine number of incorperated Rr genes and lenght of database
incCnt = 0
datCnt = 0
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        if lineLst[3] != "[]":
            incCnt += len(eval(lineLst[3]))
        datCnt += 1
    
#display results
print "\nThe results are in Smitty!:"
print "number of Rr orthos total:                              %d" % (orthCnt)
print "# of Rr orthos incorperated:                            %d" % (incCnt)
print "ratio of Rr orthos incorperated:                        %f" %\
      (float(incCnt) / float(orthCnt) * 100)
print "ratio of lines from database where Rr was incorperated: %f" %\
      (float(incCnt) / float(datCnt) * 100)


inp.close()
orth.close()
