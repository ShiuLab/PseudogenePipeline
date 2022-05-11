# This script is designed to split a karyotype file up by # symbols
# created by David E. Hufnagel on 11-04-2011

import sys

inp = open(sys.argv[1])       #input file
num = int(sys.argv[2])        #the number of files it will be split into
out1 = open(sys.argv[3], "w") #first output file
out2 = open(sys.argv[4], "w") #second output file

breakCnt = 0
#for splitting the file in 2
if num == 2:
    for line in inp:
        if breakCnt == 0:
            if line.startswith("#"):
                out1.write(line)
                breakCnt += 1
            else:
                print "***ERROR HERE***"
        elif breakCnt == 1:
            if line.startswith("#"):
                out2.write(line)
                breakCnt += 1
            else:
                out1.write(line)
        elif breakCnt == 2:
            if line.startswith("#"):
                print "***ERROR HERE***"
            else:
                out2.write(line)
        
else:
    print "***ERROR HERE***"

    
