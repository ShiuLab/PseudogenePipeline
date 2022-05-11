#This script is designed to extract At-Al, Br-Rr, single-copy-orthologs and
#At(1)-Al(1)-Br(3)-Rr(3) lines in a seperate output file for each
#Designed by David E. Hufnagel on 5-11-2012

import sys

inp = open(sys.argv[1])
atAl = open(sys.argv[2], "w")  #the output for At-Al only lines
brRr = open(sys.argv[3], "w")  #the output for Br-Rr only lines
sco = open(sys.argv[4], "w")   #the output for single-copy-ortholog lines
one23a = open(sys.argv[5], "w") #the output for At(1)-Al(1)-Br(3)-Rr(2) lines
one23b = open(sys.argv[6], "w") #the output for At(1)-Al(1)-Br(3)-Rr(1) lines

#writes the users command line prompt on the first line of the output file.
atAl.write('#python %s\n'%(' '.join(sys.argv)))
brRr.write('#python %s\n'%(' '.join(sys.argv)))
sco.write('#python %s\n'%(' '.join(sys.argv)))
one23a.write('#python %s\n'%(' '.join(sys.argv)))
one23b.write('#python %s\n'%(' '.join(sys.argv)))

#WARNING: sequence specific (to get rid of hash line)
inp.readline()

    

for line in inp:
    lineLst = line.split("\t")
    if lineLst[0] != "[]":
        At = [lineLst[0],]
    else:
        At = []
    Al = eval(lineLst[1])
    Br = eval(lineLst[2])
    Rr = eval(lineLst[3])

    #atAl extraction
    if At != [] and Al != [] and Br == [] and Rr == []:
        atAl.write(line)
    #brRr extraction
    elif At == [] and Al == [] and Br != [] and Rr != []:
        brRr.write(line)
    #sco extraction
    elif len(At) == 1 and len(Al) == 1 and len(Br) == 1 and len(Rr) == 1:
        sco.write(line)
    #one23a extraction
    elif len(At) == 1 and len(Al) == 1 and len(Br) == 3 and len(Rr) == 2:
        one23a.write(line)
    #one23b extraction
    elif len(At) == 1 and len(Al) == 1 and len(Br) == 3 and len(Rr) == 1:
        one23b.write(line)

    



inp.close()
atAl.close()
brRr.close()
sco.close()
one23a.close()
one23b.close()
