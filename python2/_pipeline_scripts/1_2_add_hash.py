#Created by: David E. Hufnagel

import sys

fd = open(sys.argv[1])
FDlist = FD.split(".")
temp = (FDlist[:-1])#.append("_hash")
temp2 = ".".join(temp)
print "temp: ", temp
print "temp2: ", temp2
temp3 = temp2 + "_hash." + FDlist[-1]
print "temp3: ", temp3

out = open(temp3,"w")

lastname = ""
cnt = 0
for line in fd:
    lineLst = line.split("\t")
    if lineLst[0] == lastname:
        out.write("\t".join(lineLst))
    elif cnt != 0:
        out.write("##\n")
        out.write("\t".join(lineLst))
    lastname = lineLst[0]
    cnt += 1
