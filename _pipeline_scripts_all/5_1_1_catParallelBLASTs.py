#This script is designed to concatenate multiple files into one file with the
#multiple files having the same names except for a number in the middle.  It
#was designed with parallel BLAST in mind.
#Created by David E. Hufnagel on Jan 4, 2012
import sys, os

prefix = sys.argv[1]    #the prefix before the number
suffix = sys.argv[2]    #the suffix after the number
high = int(sys.argv[3]) #the highest number (assumed to start at 1)



#generate line to run
outLine = "cat "
for num in range(1,high + 1):
    name = prefix + str(num) + suffix
    outLine += (name + " ")

outLine += "> %s" % (prefix.strip("_") + suffix)

#run line in unix
print outLine
os.system(outLine)
