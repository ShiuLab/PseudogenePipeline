
# int [[iL1,iR1],[iL2,iR2]..] [[pL1,pR1],[pL2,pR2],..] [e1,e2,..] ori [GI1,GI2,..]
# The file is aleady sorted based on iL1

# Modify 07/09/2009. When query overlaps, it was not dealt with.
#

import sys, math

inp = open(sys.argv[1]) # output of step2
lenT= int(sys.argv[2])
inl = inp.readline()

#a sometimes recursive function for outputing non-overlapping pseudoexons
#Ex: 2,1,1,2,1,4 --> 2 has highest e-val, output 2 --> reenter function with 4 (1 was tossed because it overlapped with 2) --> output 4
def OutputPseudogene(D):
    # get the top entry
    top = ""
    for j in D:
        if top == "" or D[j][0] > D[top][0]:
            top = j
            
    # output info
    # intergenicID, proteinGI, intergenicCoord, proteinCoord
    oup.write("%s\t%s\t%s\t%s\n" % \
          (all_contigs[i][0][0], top, str(D[top][1]), str(D[top][2])))
        
    #get the remainder (AKA the pseudoexons that were not part of or overlapping with the top e-val pseudoexon)
    topCoordLst = GetCoords(D[top][1])
    remainder = {}
    for j in D:
        currCoordLst = GetCoords(D[j][1])
        #look for overlap, if not there add the pseudoexon to remainder
        if not ((min(topCoordLst) >= min(currCoordLst) and min(topCoordLst) <= max(currCoordLst))
           or (max(topCoordLst) >= min(currCoordLst) and max(topCoordLst) <= max(currCoordLst))\
           or (min(currCoordLst) >= min(topCoordLst) and min(currCoordLst) <= max(topCoordLst))\
           or (max(currCoordLst) >= min(topCoordLst) and max(currCoordLst) <= max(topCoordLst))):
            remainder[j] = D[j]

    #if remainder not empty, send it back through this function recursively
    if list(remainder.keys()) != []:
        OutputPseudogene(remainder)
        

#get the coords for the top entry (so the min and max can be easily extracted)
def GetCoords(coords):
    coordLst = []
    for coord in coords:
        coordLst.append(coord[0])
        coordLst.append(coord[1])

    return coordLst



int1 = int2 = ""
contig  = []
all_contigs = []

print("Process pseudoexons...")
c = 0
while inl != "":
    if c % 1e3 == 0:
        print(" %i k" % (c/1e3))
    c += 1
    L = inl[:-1].split("\t")
    #print L
    int2 = L[0]
    ic2 = eval(L[1]) # intergenic seq coord
    pc2 = eval(L[2]) # protein coord
    ev2 = eval(L[3]) # evalues
    gl2 = eval(L[5]) # gene names

    #print [int1,int2,ic2,pc2,ev2,gl2]
    # new intron, process the last one
    if int1 != "":
        if int1 == int2:
            if ic2[0][0] - ic1[-1][-1] < lenT:
                contig.append([int2,ic2,pc2,ev2,gl2])
            else:
                all_contigs.append(contig)
                contig = [[int2,ic2,pc2,ev2,gl2]]
        else:
            all_contigs.append(contig)
            contig = [[int2,ic2,pc2,ev2,gl2]]
    else:
        contig = [[int2,ic2,pc2,ev2,gl2]]

    int1 = int2
    ic1 = ic2
    pc1 = pc2
    ev1 = ev2
    gl1 = gl2
    inl = inp.readline()

# Add the last one
all_contigs.append(contig)

# evaluating E value
# iterate through each pseudogene composite
print("Iterate through pseudoexon contigs...")
oup = open("%s_I%i.PS1" % (sys.argv[1],lenT),"w")
c = 0
for i in range(len(all_contigs)):
    if c % 1e3 == 0:
        print(" %i k" % (c/1e3))
    c += 1
    #print i,len(all_contigs[i])
    # D = {query_gi:[composite_e,[iC],[pC]}
    D = {}
    
    # [inter,[ic],[pc][ev][gl]]
    for j in all_contigs[i]:
        for k in range(len(j[4])):
            if j[3][k] == 0.0:  # evalue is 0
                e = 200
            else:
                e = -math.log(j[3][k],10)
            #print j[3][k],e
            if j[4][k] not in D:    # gi not in D
                D[j[4][k]] = [e,[j[1][k]],[j[2][k]]]
            else:
                D[j[4][k]][0] += e
                D[j[4][k]][1].append(j[1][k])
                D[j[4][k]][2].append(j[2][k])

    #Output info
    OutputPseudogene(D)

oup.close()
print("\npseudoexons merged!\n")
