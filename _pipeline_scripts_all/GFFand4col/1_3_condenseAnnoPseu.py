#This script is designed to condense pseudogene lines in .4col.annotated file
#with the same names. output: a .4col.annotated.condensed file

import sys

inp = open(sys.argv[1])      #input .4col.annotated file
out = open(sys.argv[2], "w") #output .4col.annotated.condensed file



#writes the users command line prompt on the first line of the output file.
out.write('#python %s\n'%(' '.join(sys.argv)))

#make dict of lines where lines with the same name are multiple values with 1 key
inp.readline()
bigDict = {}   #Dict of lines with the name as the key and the full line(s) as the value.
for line in inp:
    lineLst = line.split("\t")
    if lineLst[0] not in bigDict:
        bigDict[lineLst[0]] = [line,]
    else:
        bigDict[lineLst[0]].append(line)

#go through the dict and print new lines
for name in bigDict:
    if len(bigDict[name]) > 1:
        newline = "%s;%s" % (bigDict[name][0][:-1], bigDict[name][1].split("\t")[-1])
        out.write(newline)
    else:
        out.write(bigDict[name][0])



inp.close()
out.close()
