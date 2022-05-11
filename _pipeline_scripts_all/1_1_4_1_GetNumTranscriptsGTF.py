# This script is designed to read the string, "Cufflinks\ttranscript" to
# determine the number of transcripts in a gtf/gff file
# Created by: David E. Hufnagel 1-19-2012

import sys

inp = open(sys.argv[1]) # GFF input file

TranscriptCnt = 0
for line in inp:
    #print line
    if "Cufflinks\ttranscript" in line:
        TranscriptCnt += 1

print
print TranscriptCnt, "transcripts\n"

inp.close()
