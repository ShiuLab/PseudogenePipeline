#This script is designed to take the .comp file generated from
#R_2_B_4_linkRrBrtoAll.py and combine lines that have the same At-Al-Br. same = "123"
"""
Ex:
AT1G38440	[]	[]	[]
AT5G08160	['fgenesh2_kg.6_785_AT5G08160.1']	['Bra009325', 'Bra005973', 'Bra028670']	['RrC7653_p2_0.074']
AT5G08160	['fgenesh2_kg.6_785_AT5G08160.1']	['Bra009325', 'Bra005973', 'Bra028670']	['RrC306_p5_0.064']
AT5G08160	['fgenesh2_kg.6_785_AT5G08160.1']	['Bra009325', 'Bra005973', 'Bra028670']	['RrC69_p2_0.194']
AT4G22890	['fgenesh2_kg.7_2008_AT4G22890.3']	['Bra013653', 'Bra019349']	['RrC149_p1_0.064']

-->

AT1G38440	[]	[]	[]
AT5G08160	['fgenesh2_kg.6_785_AT5G08160.1']	['Bra009325', 'Bra005973', 'Bra028670']	['RrC7653_p2_0.074', 'RrC306_p5_0.064', 'RrC69_p2_0.194']
AT4G22890	['fgenesh2_kg.7_2008_AT4G22890.3']	['Bra013653', 'Bra019349']	['RrC149_p1_0.064']
"""
#also combines lines with the same At-Al-Rr, but different Br. same = "124"
"""
Ex:
AT1G38440	[]	[]	[]
AT5G08160	['fgenesh2_kg.6_785_AT5G08160.1']	['Bra009325']	['RrC7653_p2_0.074', 'RrC69_p2_0.194']
AT5G08160	['fgenesh2_kg.6_785_AT5G08160.1']	['Bra005973']	['RrC7653_p2_0.074', 'RrC69_p2_0.194']
AT5G08160	['fgenesh2_kg.6_785_AT5G08160.1']	['Bra028670']	['RrC7653_p2_0.074', 'RrC69_p2_0.194']
AT4G22890	['fgenesh2_kg.7_2008_AT4G22890.3']	['Bra013653', 'Bra019349']	['RrC149_p1_0.064']

-->

AT1G38440	[]	[]	[]
AT5G08160	['fgenesh2_kg.6_785_AT5G08160.1']	['Bra009325', 'Bra005973', 'Bra028670']	['RrC7653_p2_0.074', 'RrC69_p2_0.194']
AT4G22890	['fgenesh2_kg.7_2008_AT4G22890.3']	['Bra013653', 'Bra019349']	['RrC149_p1_0.064']
"""

import sys

inp = open(sys.argv[1])       #input .comp file
out = open(sys.argv[2], "w")  #temporary output .comp file
out2 = open(sys.argv[3], "w") #final output .comp file



def SaveIntoDict(gene1, gene2, dictX):
    if gene1 not in dictX:
        dictX[gene1] = [gene2]
    else:
        dictX[gene1].append(gene2)

def TheMainDeal(same):
    #writes the users command line prompt on the first line of the output file.
    out.write('#python %s\n'%(' '.join(sys.argv)))

    #WARNING: sequence specific (to get rid of hash line)
    inp.readline()

    last = []      #the last line excluding Rr
    fullLast = ""  #the last line
    radLst = []    #a temporary list of Rr sections
    for line in inp:
        lineLst = line.split("\t")
        #during same lines
        if same == "123" and lineLst[:3] == last:
            At = lineLst[0]
            Al = eval(lineLst[1])
            Br = eval(lineLst[2])
            fLastLst = fullLast.split("\t")
            lastRr = eval(fLastLst[3])
            for r in lastRr:
                radLst.append(r)
        else:
            #for the last line after same lines are done
            if radLst != []:
                At = lineLst[0]
                Al = eval(fLastLst[1])
                Br = eval(fLastLst[2])
                fLastLst = fullLast.split("\t")
                lastRr = eval(fLastLst[3])
                for R in lastRr:
                    radLst.append(R)

                radStr = str(radLst)
                newLine = "%s\t%s\t%s\t%s\n" % (At,Al,Br,radStr)
                out.write(newLine)
                
                radLst = []
            #not same lines and not after same lines
            else:
                if fullLast != "":
                    out.write(fullLast)
        if same == "123":
            last = lineLst[:3]
        fullLast = line



TheMainDeal("123")
#124 used to be here, but it was a poor approach.  I kept the format change simply because it works now, why mess with it

#the new 124 method:
inp.close()
inp = open(sys.argv[1])
#writes the users command line prompt on the first line of the output file.
out2.write('#python %s\n'%(' '.join(sys.argv)))
#WARNING: sequence specific (to get rid of hash line)
inp.readline()

bigDict = {}
for ln in inp:
    lnLst = ln.split("\t")
    sameChunk = lnLst[:2]
    temp = str(lnLst[3][:-1])
    sameChunk.append(temp)
    SaveIntoDict(str(sameChunk), str(lnLst[2]), bigDict)
    
for key in bigDict:
    k = eval(key)
    val = bigDict[key]
    #print k
    #print val
    #print val[0]
    #print
    newLine = "%s\t%s\t%s\t%s\n" % (k[0], k[1], val[0], k[2])
##    tmpLst = []
##    if k[0] == "[]" and k[1] == "[]" and len(val[0]) > 1:
##        if val[0] not in tmpLst:
##            tmpLst.append(val[0])
##        else:
##            print val[0]
    out2.write(newLine)
        

    



inp.close()
out.close()
out2.close()
