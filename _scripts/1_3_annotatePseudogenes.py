#This script is designed to extract pseudogenes from a .disable_count file,
#anjust the coordinates so they are genomic coordinates, and export that info
#into a 4col output file.
#Designed by David E. Hufnagel on May 24, 2012
"""Algorithm:
1) Go through pseudogene file, extract and alter coordinates, extract
ID(whole original line joined by ";"s) and four # identifier.
2) Write the output into a new 4 column file"""
import sys

pseu = open(sys.argv[1])     #input disable_count file with pseudogene info
out = open(sys.argv[2], "w") #output .4col file



out.write('#python %s\n'%(' '.join(sys.argv))) #writes the users command line prompt on the first line of the output file.

for line in pseu:
    lineLst = line.split(" ")
    if not line.startswith("\n"):
        if line.startswith("#"):
            oldPsCoords = tuple(lineLst[2].split(":")[1].split("-"))
            oldGenCoords = tuple(lineLst[1].split("|")[-1].split("-"))
            newCoords = ((int(oldGenCoords[0]) + int(oldPsCoords[0])) - 1, \
                         (int(oldGenCoords[0]) + int(oldPsCoords[1])) - 1)
            name = (";".join(lineLst[:3]) + ";" + "|".join(lineLst[4:])[:-1])[1:]
            chro = "|".join(lineLst[1].split("|")[:-1])
            newline = "%s\t%s\t%s\t%s\n" % (name, chro, newCoords[0], newCoords[1])
            out.write(newline)
        





pseu.close()
out.close()
