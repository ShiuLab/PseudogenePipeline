#This script is designed to take a 3 column gff file with the first column
#containing "species|chromosome" and making several 4 column gff files in the
#format "annotation chromosome startCoord stopCoord"
#Created by David E. Hufnagel on June 4, 2012
#WARNING HIGHLY SITUATION SPECIFIC

import sys

inp = open(sys.argv[1])
one = open("ChlNC64A_1.4colSimple", "w")
two = open("Coc_C169_1.4colSimple", "w")
three = open("Cre.4colSimple", "w")
four = open("MicpuC2.4colSimple", "w")
five = open("MicpuN3.4colSimple", "w")
six = open("Ost9901_3.4colSimple", "w")
seven = open("OstRCC809_2.4colSimple", "w")
eight = open("Ostta4.4colSimple", "w")
nine = open("Vcarteri.4colSimple", "w")



for line in inp:
    lineLst = line.split("\t")
    newLine = "annotation\t%s\t%s\t%s" % (lineLst[0].split("|")[1],lineLst[1], lineLst[2])
    if lineLst[0].split("|")[0] == "ChlNC64A_1":
        one.write(newLine)
    elif lineLst[0].split("|")[0] == "Coc_C169_1":
        two.write(newLine)
    elif lineLst[0].split("|")[0] == "Cre":
        three.write(newLine)
    elif lineLst[0].split("|")[0] == "MicpuC2":
        four.write(newLine)
    elif lineLst[0].split("|")[0] == "MicpuN3":
        five.write(newLine)
    elif lineLst[0].split("|")[0] == "Ost9901_3":
        six.write(newLine)
    elif lineLst[0].split("|")[0] == "OstRCC809_2":
        seven.write(newLine)
    elif lineLst[0].split("|")[0] == "Ostta4":
        eight.write(newLine)
    elif lineLst[0].split("|")[0] == "Vcarteri":
        nine.write(newLine)
    else:
        print "ERROR LINE: ", lineLst[0].split("|")[0]
    







inp.close()
one.close()
two.close()
three.close()
four.close()
five.close()
six.close()
seven.close()
eight.close()
nine.close()
