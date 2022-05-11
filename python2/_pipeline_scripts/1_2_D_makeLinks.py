#This script is designed to make a CIRCOS links file from a specific est file:
"""
#EST	LinkageGrp		CentiMorg	Radish	Brassica	Lyrata	Thaliana
RSCR04O08	LG1	0	0	RrC20787:1292-2009	A02:789448-790172	7:21246447-21247083	Chr5:16738940-16739513
RSCS14I08	LG1	16000000	16	RrC9941:1452-2161	A02:1007920-1008694	5:16439043-16439496	Chr3:19688308-19688734
-->
AthAly000	Chr1	930828	3255693
AthAly000	scaffold_1	1125712	3818257

Link001	Chr1	930828	3255693
Link001	scaffold_1	1125712	3818257"""
#Created by David E. Hufnagel on May 22, 2012
#WARNING: HIGHLY SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])      #input est file
tmp = open("temp.temp", "w")
out = open(sys.argv[2], "w") #output .links file





def FixChromo(spe, chromo):
    if spe == "Br":
        newChromo = "Br" + chromo
    elif spe == "Al":
        newChromo = "Al" + chromo
    elif spe == "At":
        newChromo = "At" + chromo

    return newChromo
    
def ProcessChunk(spe, RrLine, chunk):
    chromo = chunk.split(":")[0]
    bad = False
    newLines = ""
    
    if spe == "Br":
        if chromo.startswith("S"):
            bad = True
        chromo = FixChromo("Br", chromo)
        coordEnd = chunk.split(":")[1].split("-")[1]
    elif spe == "Al":
        chromo = FixChromo("Al", chromo)
        coordEnd = chunk.split(":")[1].split("-")[1]
    elif spe == "At":
        chromo = FixChromo("At", chromo)
        coordEnd = chunk.split(":")[1].split("-")[1][:-1]
        
    coordStart = chunk.split(":")[1].split("-")[0]

    if bad == False:
        newLines = "%s%s\t%s\t%s\n" % (RrLine, chromo,\
                                        coordStart, coordEnd)
    return newLines

def AddZeroes(num): # 3  -> "003"
    tempNum = str(num)
    while len(tempNum)<3:
        tempNum = "0" + tempNum    
    
    newNum = tempNum
    return newNum







inp.readline()
for line in inp:
    lineLst = line.split("\t")
    RrLine = "Rr%s\t%s\t%d\n" % (lineLst[1], lineLst[2],\
                                   (int(lineLst[2]) + 1000))

    BrName = lineLst[5]  #just the chromosome name and coords
    newLinesBr = ProcessChunk("Br", RrLine, BrName)

    AlName = lineLst[6]
    newLinesAl = ProcessChunk("Al", RrLine, AlName)

    AtName = lineLst[7]
    newLinesAt = ProcessChunk("At", RrLine, AtName)

    if newLinesBr != "":
        tmp.write(newLinesBr)
        tmp.write(newLinesAl)
        tmp.write(newLinesAt)
        

#adds linknums to each line
tmp.close()
tmp = open("temp.temp")


#writes the users command line prompt on the first line of the output file.
#out.write("#python %s\n"%(" ".join(sys.argv)))

cnt = 0
cnt2 = 1
for line in tmp:
    linkNum = "link" + AddZeroes(cnt2)
    newLine = "%s\t%s" % (linkNum, line)
    out.write(newLine)
    cnt += 1

    if not cnt % 2:
        cnt2 += 1



inp.close()
tmp.close()
out.close()
