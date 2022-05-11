#take a phylip .mlout file and remove all info, but the trees
#Created By David E. Hufnagel on Oct 11, 2012

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



doWrite = False
for line in inp:
    if line.startswith("(Transition/transversion parameter ="):
        doWrite = True
    elif line.startswith("remember:"):
        doWrite = False
        
    if doWrite == True:
        if not line.startswith("(Transition/transversion parameter ="): #don't print start delimiter line
            out.write(line)
        else:
            out.write("####")





inp.close()
out.close()
