# This script is intended to filter a GMAP output file to keep only
# results passing certain thresholds.
# Created by David E. Hufnagel on July 30, 2012 (based on 1_2_4_3_FilterBlast.py)

import sys

fd = open(sys.argv[1])        # The GMAP output filename
iden = float(sys.argv[2])     # The percent ID cutoff (will keep this value or bigger)
mlen = int(sys.argv[3])       # The minimum length cutoff
out = open(sys.argv[4], "w")  # The output filename





#print the unix command line that called this script
out.write('#python %s\n'%(' '.join(sys.argv)))

for line in fd:
    if not line.startswith("#"):
        doWrite = True
        if not line.startswith("\n"):
            lineLst = line.split("\t")
            
            # % ID filter
            if (not float(lineLst[5]) >= iden) and (iden != 0.0):
                doWrite = False
            # Min length filter
            if (not int(lineLst[2]) >= mlen) and (mlen != 0):
                doWrite = False

            if doWrite == True:
                out.write(line)

fd.close()
out.close()
