#This script is designed to make a number based reference file (starting with 0)
#from names in a fasta file.  It was originally designed for handling Vmatch out
#put which gives a sequence number instead of name (morons!)
#Created by David E. Hufnagel on July 17, 2012

import sys

inp = open(sys.argv[1])      #Input tabular file with the names to make the reference file from
ref = open(sys.argv[2], "w") #Output reference file with numbers in the first column and names in the second





#Write the users command line prompt on the first line of the output file.
ref.write("#python %s\n" % (" ".join(sys.argv)))

num = 0
for line in inp:
    if not line.startswith("#"):
        if line.startswith(">"):
            name = line[1:].strip()
            newLine = "%s\t%s\n" % (num, name)
            ref.write(newLine)
            num += 1






inp.close()
ref.close()
