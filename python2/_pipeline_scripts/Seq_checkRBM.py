#This script is designed to double check that R_2_A_1_getRecipicalBlast.py
#did it's job correctly by checking the .rbm file for the total number of pairs
#and the number of each species in each column
#WARNING: ONLY HAS CAPACITY FOR At, Al, Br AND Rr THUS FAR
"""theory:
At = AT
Al = fg, sc, Al
Rr = Rr
Br = Br"""
import sys

inp = open(sys.argv[1])      #the input .rbm file to be checked
out = open(sys.argv[2], "w") #the output file with the data results



#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

pairCnt = 0
AtCol1Cnt = 0
AlCol1Cnt = 0
BrCol1Cnt = 0
RrCol1Cnt = 0
AtCol2Cnt = 0
AlCol2Cnt = 0
BrCol2Cnt = 0
RrCol2Cnt = 0
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        
        if lineLst[0].startswith("AT"):
            AtCol1Cnt += 1
        elif lineLst[0].startswith("fg") or lineLst[0].startswith("sc") or lineLst[0].startswith("Al"):
            #print lineLst[0]
            AlCol1Cnt += 1
        elif lineLst[0].startswith("Br"):
            #print lineLst[0]
            BrCol1Cnt += 1
        elif lineLst[0].startswith("Rr"):
            RrCol1Cnt += 1
        else:
            print "error1: ", line

        if lineLst[1].startswith("AT"):
            AtCol2Cnt += 1
        elif lineLst[1].startswith("fg") or lineLst[1].startswith("sc") or lineLst[1].startswith("Al"):
            #print lineLst[1]
            AlCol2Cnt += 1
        elif lineLst[1].startswith("Br"):
            #print lineLst[1]
            BrCol2Cnt += 1
        elif lineLst[1].startswith("Rr"):
            RrCol2Cnt += 1
        else:
            print "error2: ", line
        
        pairCnt += 1

out.write("\ntotal: %s\n" % (pairCnt))
out.write("At1: %s\n" % (AtCol1Cnt))
out.write("Al1: %s\n" % (AlCol1Cnt))
out.write("Br1: %s\n" % (BrCol1Cnt))
out.write("Rr1: %s\n" % (RrCol1Cnt))
out.write("At2: %s\n" % (AtCol2Cnt))
out.write("Al2: %s\n" % (AlCol2Cnt))
out.write("Br2: %s\n" % (BrCol2Cnt))
out.write("Rr2: %s\n" % (RrCol2Cnt))
    
    





inp.close()
out.close()
