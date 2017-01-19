# This script is intended to simply filter a blast output file to keep only
# results passing certain thresholds.
# Created by David E. Hufnagel
# Updated 01-04-2012: added filters for % coverage and minimum length
#   (only tested for 12 column files).  Any unwanted variable should be entered
#   as 0.  Variables will be printed off
# note: for percent coverage you will need to add a .size file that is simply
# a list of the gene names in the query with the length of the seqs (this can
# most easily be acquired from the FastaManager.py function, get_size)
#   Ex: CUFF.1\t72\nCUFF.2\t312...
# note: coverage will only work for blastn or blastp without altering this script

import sys

fd = open(sys.argv[1])        # The blast filename
iden = float(sys.argv[2])     # The percent ID cutoff (will keep this value or bigger)
e = float(sys.argv[3])        # The e-value cutoff (will keep this value or smaller)
mlen = int(sys.argv[4])       # The minimum length cutoff
col = int(sys.argv[5])        # The number of columns in the file (handles 3 or 12 (m8))
out = open(sys.argv[6], "w")  # The output filename





#print off the stated variables so the user can double check it
print "\nyour variables:"
print "input file = %s" % (fd)
if iden != 0.0:
    print "minimum percent ID = %s" % (iden)
if e != 0.0:
    print "maximum e-value = %s" % (e)
if mlen != 0:
    print "minimum length = %s" % (mlen)
if col != 0:
    print "number of blast columns = %s" % (col)
print "output file = %s\n" % (out)

#print the unix command line that called this script
out.write('#python %s\n'%(' '.join(sys.argv)))

#WARNING: NOT RECENTLY TESTED
# For 3 column blast files (gene1   gene2   e-value) 
if col == 3:
    for line in fd:
        lineLst = line.split("\t")
        if float(lineLst[2]) <= e:
            out.write(line)

# For 12 column blast files (Query id   Subject id   % identity alignment
#length   mismatches   gap   openings   q. start   q. end s. start   s. end
#e-value   bit score)
elif col == 12:
    for line in fd:
        if not line.startswith("#"):
            doWrite = True
            if not line.startswith("\n"):
                lineLst = line.split("\t")
                
                # % ID filter
                if (not float(lineLst[2]) >= iden) and (iden != 0.0):
                    doWrite = False
                # E-value filter
                if (not float(lineLst[10]) <= e) and (e != 0.0):
                    doWrite = False
                # Min length filter
                if (not int(lineLst[3]) >= mlen) and (mlen != 0):
                    doWrite = False

                if doWrite == True:
                    out.write(line)

fd.close()
out.close()
