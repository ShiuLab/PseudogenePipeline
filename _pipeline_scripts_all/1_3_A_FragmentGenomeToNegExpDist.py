#This script is designed to fragment a genome (At) in a fashion that creates a
#negative exponential distribution with certain parameters (like the Rr genome)
#Created by David E. Hufnagel on July 25, 2012
#WARNING: INCOMPLETE
"""Algorithm:"""


import sys, random

mean = 1.0 / float(sys.argv[1])
num = int(sys.argv[2])
out = open(sys.argv[3], "w")



for x in range(1, num+1):
    y = random.expovariate(mean)
    newLine = str(y) + "\n"
    out.write(newLine)




out.close()
