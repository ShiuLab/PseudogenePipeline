#This script was designed to make a 4col file from a gff file.  This is a
#dynamic script and will likely be changed many times to fit particular
#situations.
#Created by David E. Hufnagel on Dec 16, 2012

import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")



for line in inp:
    if not line.startswith("#"):
        lineLst = line.strip().split("\t")
        if lineLst[2] == "gene":
            name = lineLst[8].split("Name=")[1].strip(";")[:]
            #name = "Vcarteri|" + name.strip(";")
            #if name.startswith("GR"):
            #    name = "".join(name.split("_")[:-1])
            #pyname = name.replace("PGSC0003DMG","PGSC0003DMP")
            newLine = "%s\t%s\t%s\t%s\n" % (name, lineLst[0], lineLst[3], lineLst[4])
            out.write(newLine)


inp.close()
out.close()
