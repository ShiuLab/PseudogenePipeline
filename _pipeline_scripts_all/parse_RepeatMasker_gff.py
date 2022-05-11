import sys
print "INP1: RepeatMasker *.out output"
print "INP2: Cutoff SW score"
print "INP3: Cutoff Divergence (20.0,15.0,10.0...)"
file1=open(sys.argv[1],'r')
thresh=int(sys.argv[2])
div=float(sys.argv[3])
out1=open(sys.argv[1]+".%iCutoff%.1fdivg.4col"%(thresh,div),'w')

line1=file1.readline()
m=0;d=0;dict1={}
while line1:
    tab1=line1.strip().split()
    if tab1!=[]:
        i1=tab1[0]
        if i1.startswith('#') or i1.startswith('SW') or i1.startswith('score'):
            pass
        else:
            #print tab1
            d+=1
            score=int(tab1[0]); divergence=float(tab1[1])
            if score>=thresh and divergence<=div:
                m+=1                
                chr1=tab1[4]; st=int(tab1[5]); end=int(tab1[6]); o=tab1[8]
                                
                #print chr1, st, end, o
                mtp1=tab1[10];mtp2=tab1[9]
                motif=('%s|%s'%(mtp1,mtp2))
                if o=='+':
                    v1=st;v2=end
                elif o=='-':
                    v1=end;v2=st
                else:
                    v1=st;v2=end
                name=('%s|%s-%s|%s|%s'%(chr1,v1,v2,score,divergence))
                out1.write('%s\t%s\t%s\t%s\t%s\n'%(name,chr1,v1,v2,motif))

                #Parse Motif
                mtype1=motif.split('|')[0]
                if '/' in mtype1:
                    topType=mtype1.split('/')[0]
                    secondType=mtype1.split('/')[1]
                else:
                    topType=mtype1
                    secondType=mtype1
                if topType not in dict1:
                    dict1[topType]={}
                    dict1[topType][secondType]=[motif]
                else:
                    if secondType not in dict1[topType]:
                        dict1[topType][secondType]=[motif]
                    else:
                        dict1[topType][secondType].append(motif)                        
    line1=file1.readline()
file1.close()
out1.close()

out2=open(sys.argv[1]+".%iCutoff%.1fdivg.types"%(thresh,div),'w')
for v1 in dict1:
    for v2 in dict1[v1]:
        count=len(dict1[v1][v2])
        out2.write('%s\t%s\t%s\n'%(v1,v2,count))
        
out2.close()        
print "Total features: ", d
print "No. of features written to 4col file: ", m
print "Done!"
    
