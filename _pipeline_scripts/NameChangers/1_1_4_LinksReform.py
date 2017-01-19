
import sys

in1 = open(sys.argv[1])      #input filename
out = open(sys.argv[2], "w") #output filename



#constants
N1Start = "Chr"
N2Start = "scaffold_"
N1End = "At"
N2End = "Al_scaf_"
ScafNumLen = 4



def AddZeroes(num): # 3  -> "0003"
    tempNum = str(num)
    while len(tempNum)<ScafNumLen:
        tempNum = "0" + tempNum    
    
    newNum = tempNum
    return newNum



for line in in1:
    print line
    lineLst = line.split("\t")
    if lineLst[1].startswith(N1Start):
        newName = N1End + lineLst[1][3:]
    elif lineLst[1].startswith(N2Start):
        temp = lineLst[1].split("_")
        temp2 = AddZeroes(temp[1])
        newName = N2End + temp2
    else:
        print "\n\n\n\n***ERROR***\n\n\n\n"
    print

    newLine = "%s\t%s\t%s\t%s\n" % (lineLst[0], newName, lineLst[2], lineLst[3].strip())
    out.write(newLine)

in1.close()
out.close()
