#This script was designed to take an unknown number of input cuffdiff files and
#making and edgeR input table from the information in the format:
#geneName   file1ReadCount   file2ReadCount   file3ReadCount
#Created by David E. Hufnagel on Oct 29, 2012
#WARNING: THIS SCRIPT ASSUMES ALL INPUT FILES HAVE THE SAME GENE LIST
"""Algorithm:
1) Go through input files.  For each input file make a dict of genes for each
   file with key: name val: count
2) output info by going through one dict and pulling in info from all the
   others"""

import sys

cufDifs = sys.argv[1:-1]       #the input cuffdif output files
out = open(sys.argv[-1], "w")  #output edgeR input file




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))
#output title
out.write("#")
for fileName in cufDifs:
    out.write("  %s" % (fileName))
out.write("\n")
    
#1) Go through input files.  For each input file make a dict of genes for each
#   file with key: name val: count
dictLst = []
for fileName in cufDifs:
    fd = open(fileName)
    fileDict = {}
    for line in fd:
        lineLst = line.strip("\n").split("\t")
        fileDict[lineLst[0]] = lineLst[1]

    dictLst.append(fileDict)

#2) output info by going through one dict and pulling in info from all the others
for name, cnt in dictLst[0].items():
    out.write("%s\t%s" % (name, cnt))
    for dicty in dictLst[1:]:
        out.write("\t%s" % (dicty[name]))
    out.write("\n")
            




out.close()
