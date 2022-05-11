#This script was designed to convert a .PaCE file to a .pillars file
"""
{Cluster#}   47898

  {Member#}  >gi|167504305|gb|FD989781.1|FD989781

  {Member#}  >gi|154091389|gb|EV524362.1|EV524362

-->

gi|167504305|gb|FD989781.1|FD989781     gi|154091389|gb|EV524362.1|EV524362
"""

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")





#Go through inp and extract clusters and names into a list of lists
bigLst = []
tmpLst = []
for line in inp:
    #print line
    if line.startswith("{Cluster#}"):
        if tmpLst != []:
            bigLst.append(tmpLst)
        tmpLst = []
    elif "{Member#}" in line:
        name = line.split(">")[1].strip()
        name = name.replace(".","_")
        tmpLst.append(name)
else:
    bigLst.append(tmpLst)
        
#Go through bigLst and write newLines
for cluster in bigLst:
    newLine = "\t".join(cluster) + "\n"
    out.write(newLine)









        



inp.close()
out.close()
