
import sys,re

def space(s):
    L = s.split(" ")
    tmp = []
    for i in L:
        if i != "":
            tmp.append(i)
    return tmp

def finditer(p,s):
    p = re.compile(p)   # pattern with one or more "-"
    m = re.finditer(p,s) # iterative find all p
    spans = []
    for i in m:
        spans.append(i.span())
    return spans
    
# From: http://www.d.umn.edu/~mhampton/m5233s7.html
def nw(seq1, seq2, matrix, gap = -8, match = 1, mismatch = -1):
    
    # initialize scoring and 'arrow' matrices to 0
    S = [[0 for x in range(len(seq2)+1)] for y in range(len(seq1)+1)]
    P = [[0 for x in range(len(seq2)+1)] for y in range(len(seq1)+1)]
    
    # initialize borders
    # for P (arrows), use 2 for diagonal, -1 for horizontal, and 1 for 
    # vertical moves.
    # I have tried to consistently use i for rows (vertical positions) in the 
    # score and pointer tables, and j for columns (horizontal positions).
    for i in range(len(seq1)+1):
        S[i][0] = gap*i
        P[i][0] = 1
    for j in range(len(seq2)+1):
        S[0][j] = gap*j
        P[0][j] = -1
    # fill with scores
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            L1 = seq1[i-1]
            L2 = seq2[j-1]
            # treat stop or frameshift as gap
            # | is the substitute for \ because of it is also use as escape
            if L1 in ["*","/","\\","|"]:
                L1 = "-"
            if L2 in ["*","/","\\","|"]:
                L2 = "-"
            Sdg = S[i-1][j-1] + matrix[L1][L2]  # diagnoal
            Shz = S[i][j-1] + gap              # horizontal
            Sup = S[i-1][j] + gap              # up
            # TempS is list of the three S and their P
            TempS = [[Sdg,2],[Shz,-1],[Sup,1]]
            # Now we keep the highest score, and the associated direction (pointer)
            S[i][j], P[i][j] = max(TempS)
    # backtrace from end
    [i,j] = [len(seq1),len(seq2)]
    align1 = ''
    align2 = ''
    # does NOT traceback on multiple optimal paths.
    while [i,j] != [0,0]:
        if P[i][j] == 2:
            align1 = align1 + seq1[i-1]
            align2 = align2 + seq2[j-1]
            i = i - 1
            j = j - 1
        elif P[i][j] == -1:
            align1 = align1 + '-'
            align2 = align2 + seq2[j-1]
            j = j - 1
        else:
            align1 = align1 + seq1[i-1]
            align2 = align2 + '-'
            i = i - 1
    # the alignments have been created backwards, so we need to reverse them:
    align1 = align1[::-1]
    align2 = align2[::-1]
    # print out alignment
    return align1,align2

print("Read BLOSUM50 matrix...")
# read blosum50 matrix
inp = open(sys.argv[1])
inl = inp.readlines()
aa  = space(inl[0][:-1])
b50 = {}
inl = inl[1:]
for i in range(len(inl)):
    #print [inl[i]]
    s = space(inl[i][:-1])
    for j in range(len(aa)):
        #print aa[j],
        if aa[j] not in b50:
            b50[aa[j]] = {aa[i]:int(s[j])}
            #print "new:",b50[aa[j]]
        else:
            b50[aa[j]][aa[i]] = int(s[j])
            #print "app:",b50[aa[j]]

print("Read the sw.out file...")
# As of 2/16,2012, change the following so it will only allow one alignment
# to be reported
inp = open(sys.argv[2])
inl = inp.readline()
pair = []
stat = []
seq  = []
sw   = {} # {pair:[[stat1,seq1,stop_1,stop2,frame1,frame2],...]}
          # stop1 and frame1: no alignment necessary
          # stop2 and frame2: aligned because they are too close to gaps

while inl != "":
    if inl[0] == "#":
        pair = inl[1:-1]
        stat = []
    elif inl[0] == ">":
        stat = inl[1:-1].split(" ")
        seq  = []
    else:
        # rid of the escape character and replace it with "2"
        L = inl[:-1].split("\\")
        L = "|".join(L)
        seq.append(L)
    inl = inp.readline()

    # Modified to allow only one entry, 2/16/2012
    if stat != [] and len(seq) == 2:
        if pair not in sw:
            sw[pair] = [[stat,seq,0,0,0,0]]
        # to prevent redundant entries, 8/12, 07'
        #elif [stat,seq,0,0,0,0] not in sw[pair]:
        #   sw[pair].append([stat,seq,0,0,0,0])
    
# Now compare the sequences
print("Compare sequences:")
p1   = "[-]+"  # any number of "-"
p2   = "[*]"   # a single stop
p3   = "[/|]"  # single "/" or "\"
flank= 5       # need stop and frameshift to be be this number of amino acids
               # away from gap with length > 1
gaplen = 1   # need gap to be larger than this to be considered
oup = open(sys.argv[2]+".disable_count","w")
c = 0
swkl = len(list(sw.keys()))
print(" total: %i alignments" % swkl)
for i in sw:
    if c % 1e3 == 0:
        print("",c/1e3)
    c += 1
    #print "",i,len(sw[i])
    c = 0
    for j in sw[i]:
        # get gap list in seq1
        gap1 = finditer(p1,j[1][0])
        for k in range(len(gap1)):
            gap1[k] = [gap1[k][0],gap1[k][1],"-"]
        #print gap1
        
        # get stop list in seq2
        stop2 = finditer(p2,j[1][1])
        for k in range(len(stop2)):
            stop2[k] = [stop2[k][0],stop2[k][1],"*"]
        #print stop2
        
        # get framshift list in seq2
        fram2 = finditer(p3,j[1][1])
        for k in range(len(fram2)):
            fram2[k] = [fram2[k][0],fram2[k][1],"f"]        
        #print fram2
        
        c1 = stop2 + fram2
        c1.sort()
        #print c1
        # combine check stop and frameshift
        for k in c1:
            #print k
            stop_ok1 = 1
            stop_ok2 = 0
            fram_ok1 = 1
            fram_ok2 = 0
            for m in gap1:
                # gap has to be larger than the threshold
                if m[1]-m[0] > gaplen:
                    # x-----xxxx
                    # xxxx*xxxxx
                    #print k,m
                    if k[0] >= m[0] and k[1] < m[1]:
                        if k[2] == "*":
                            stop_ok1 = 0
                        else:
                            fram_ok1 = 0
                        break
                    #     12345
                    # x----xxxxxx
                    # xxxxxxx*xxx
                    # minus 1 because gap end is exclusive
                    elif k[0] >= m[0] and k[0] < m[1]-1+flank:
                        # do ungapped global alignment
                        seq1 = j[1][0][m[1]:k[1]+2]
                        seq2 = j[1][1][m[0]:k[1]+2]
                        #print "1:",seq1,seq2
                        aln1,aln2 = nw(seq1,seq2,b50)
                        #print aln1
                        #print aln2
                        # evaluate if the stop/frameshift is sitting in a gap
                        # stop/frameshift position in the alignment: k[0]-m[1]
                        sf = k[0]-m[0]
                        if aln1[sf] == "-" and \
                           (aln1[sf-1] == "-" or aln1[sf+1] == "-"):
                            if k[2] == "*":
                                stop_ok1 = 0
                            else:
                                fram_ok1 = 0
                            break
                        else:
                            if aln2[sf] == "*":
                                stop_ok2 = 1
                            else:
                                fram_ok2 = 1
                        
                        #print aln2[sf],stop_ok2,fram_ok2
                    #  12345
                    # xxxxxx---xx
                    # xxx*xxxxxxx
                    elif k[0] < m[0] and k[0]+flank >= m[0]:
                        # do ungapped alignment
                        # make sure left coord >= 0
                        sL = max(k[0]-2,0)
                        seq1 = j[1][0][sL:m[0]]
                        seq2 = j[1][1][sL:m[1]]
                        #print "2:",seq1,seq2
                        aln1,aln2 = nw(seq1,seq2,b50)
                        #print aln1
                        #print aln2
                        sf = k[0]
                        if k[0]-2 >= 0:
                            sf = 2
                            
                        if aln1[sf] == "-" and \
                           (aln1[sf-1] == "-" or aln1[sf+1] == "-"):
                            if k[2] == "*":
                                stop_ok1 = 0
                            else:
                                fram_ok1 = 0
                            break
                        else:
                            if aln2[sf] == "*":
                                stop_ok2 = 1
                            else:
                                fram_ok2 = 1
            if k[2] == "*" and stop_ok1:
                if stop_ok2 == 0:
                    #print "add1"
                    sw[i][c][2] += 1 # stop that are relatively far from gap
                else:
                    #print "add2"
                    sw[i][c][3] += 1 # stop qualified after realignment
            elif k[2] != "*" and fram_ok1:
                if fram_ok2 == 0:
                    #print "add3"
                    sw[i][c][4] += 1 # frame shift that are far from gap
                else:
                    #print "add4"
                    sw[i][c][5] += 1 # frame shift qualified after realignment
        
        
        # generate output
        # #plantp interg range STOP1=X STOP2=X FS1=X FS2=X
        # seq1
        # seq2
        oup.write("#%s %s %s %i %i %i %i\n" % \
                        (i,sw[i][c][0][-1],sw[i][c][0][3],sw[i][c][2],
                         sw[i][c][3],sw[i][c][4],sw[i][c][5]))
        oup.write("%s\n%s\n\n" % (sw[i][c][1][0],sw[i][c][1][1]))
    
        #print sw[i][c]
        c += 1

oup.close()
print("Done!")
