#This script is designed to create a file for the Fisher Exact test in the
#format name\t1\t2\t3\t4 where 1=spe1+ 2=spe1- 3=spe2+ 4=spe2- from a GO slim file
#Created by David E. Hufnagel on July 12, 2012

import sys

inp1 = open(sys.argv[1])    #input GOslim file for spe1
inp2 = open(sys.argv[2])    #input GOslim file for spe2
out = open(sys.argv[3], "w") #output fisher-exact file





def ProcessInput(fd):
    inpDict = {}
    inpTot = 0
    for line in fd:
        if not line.startswith("#"):
            lineLst = line.split("\t")
            name = "_".join(lineLst[1].split(" "))
            inpDict[name] = lineLst[2].strip()
            inpTot += int(lineLst[2])

    return inpDict, inpTot



        
    
#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through inp1 and make a dict (key: GO name val: num) and get a total
#of nums from inp1
inp1Dict, inp1Tot = ProcessInput(inp1)
#Go through inp2 and make a dict of key: GO names val: num and get a total
#of nums from inp2
inp2Dict, inp2Tot = ProcessInput(inp2)

#Go through the 2 dicts and output the information
for key1 in inp1Dict:
    if key1 in inp2Dict:
        val1 = int(inp1Dict[key1])
        val2 = int(inp2Dict[key1])
        newLine = "%s\t%s\t%s\t%s\t%s\n" % (key1, val1, (inp1Tot - val1), val2, (inp2Tot - val2))
        out.write(newLine)
    else:
        val1 = int(inp1Dict[key1])
        val2 = 0
        newLine = "%s\t%s\t%s\t%s\t%s\n" % (key1, val1, (inp1Tot - val1), val2, (inp2Tot - val2))
        out.write(newLine)


inp1.close()
inp2.close()
out.close()
