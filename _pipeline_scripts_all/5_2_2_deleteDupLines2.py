#This script is designed to delete duplicate lines in a file when
#they're next to each other
#Created by David E. Hufnagel on April 4, 2013
import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")

lastLine = ""
cnt = 0
for line in inp:
    if line not == lastLine:
        lastLine = line
        out.write(line)

    if not cnt % 1000:
        print cnt, " lines processed"
    cnt += 1


inp.close()
out.close()
