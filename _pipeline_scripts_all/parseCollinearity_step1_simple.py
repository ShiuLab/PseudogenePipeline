import sys
print "INP1: MCSCANX collinearity output"

file1=open(sys.argv[1],'r')
out1=open(sys.argv[1]+".tab",'w')
out1.write('#python %s\n'%(' '.join(sys.argv)))
out1.write('#G1\tG2\tBlockID|PairID|EVALUE\n')
line1=file1.readline()
dict1={}; dict2={}; m=0
dict1['AT']=1; dict1['AL']=2; dict1['Br']=3; dict1['Rr']=4
dict1['Ps']=5

while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        add1=0; add2=0
        if tab1!=['']:
            block=tab1[0].split('-')[0]
            pair=tab1[0].split('-')[1].split()[0].split(':')[0]
            ev=tab1[-1].split()[0]            
            val3=('%s|%s|%s'%(block,pair,ev))
            g1=tab1[1]; g2=tab1[2];
            id1=g1[0:2]; id2=g2[0:2]
            
            out1.write('%s\t%s\t%s\n'%(ng1,ng2,val3))            
            add1=0; add2=0

            if nid1 not in dict2:
                dict2[nid1]={}
                dict2[nid1][nid2]=1
            else:
                if nid2 not in dict2[nid1]:
                    dict2[nid1][nid2]=1
                else:
                    dict2[nid1][nid2]+=1
        m+=1        
    line1=file1.readline()
file1.close(); out1.close()

for id1 in dict2:
    for id2 in dict2[id1]:
        print id1, id2, dict2[id1][id2]



print "Done!"

