#This script was designed to make a 4col file from a gff file.  This is a
#dynamic script and will likely be changed many times to fit particular
#situations.  This second version is designed to handle gff1 format
#Created by David E. Hufnagel on May 29, 2013

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "gene":
            #name = lineLst[8].split("name")[1].split(";")[0].strip().strip('"')
            name= lineLst[8].split("transcriptId")[1].split(";")[0].strip()
            name = "Ostta4|" + name
            newLine = "%s\t%s\t%s\t%s\n" % (name, lineLst[0], lineLst[3], lineLst[4])
            out.write(newLine)


inp.close()
out.close()
