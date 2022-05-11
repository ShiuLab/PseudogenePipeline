#This script is designed to do all things related to changing names within a
#file.  Example: in col 1 AT3G34029.1 --> AT3G34029
#Created by David E. Hufnagel on July 9, 2012

import sys

class NameChanger:
    def Help(self):
        print "\nParameters:"
        print "    file1 - the first file for comparison."
        print "    file2 - the second file for comparison."
        print "    o - ordered.  Whether the two files should be compared"
        print "        line-by-line in order, y or n"
        print "    p - print.  Whether the differing lines should be printed, y or n"



NC = NameChanger()

for i in range(1,len(sys.argv),2):
    if sys.argv[i] == "-input":
        file1 = open(sys.argv[i+1])
    elif sys.argv[i] == "-output":
        file2 = open(sys.argv[i+1])
    elif sys.argv[i] == "-format":
        o = sys.argv[i+1]
    elif sys.argv[i] == "-spe":
        p = sys.argv[i+1]
    elif sys.argv[i] == "-h":
        NC.Help()
    else:
        print "UNKNOWN FLAG:",sys.argv[i]
        print "add -h to get help."
        sys.exit(0)

##if o == "y":
##    Comp.Ordered(num)
##elif o == "n":
##    Comp.Unordered(num)
##else:
##    print "improper -o.  Should be 'y' or 'n'"
##    print "add -h to get help."
        
