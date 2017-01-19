import sys
print "INP1: Seq file of all AT-AL sequences"
print "INP2: Pillars file with ortholog IDs"
file1=open(sys.argv[1], 'r').readlines()
dict1={}
m=0
n=0
for i in range(0,len(file1)-1):
    line1=file1[i]
    if line1.startswith('>')==True:
        #g1=line1.strip()[1:].split('|')[0].split('.')[0]
        g1=line1.strip()[1:]
        seq=file1[i+1]        
        if g1 not in dict1:
            dict1[g1]=seq

file2=open(sys.argv[2], 'r')
out1=open(sys.argv[2]+".fa", 'w')
out1.write('#python %s\n'%(' '.join(sys.argv)))
line2=file2.readline()
dict2={}; x=0; f=0
while line2:
    if line2.startswith('#')==True:
        pass
    else:
        tab2=line2.strip().split('\t')
        tlen=len(tab2)
        i=0
        for item in tab2:
            #nlist=[]
            #nlist.append(item)            
            #nitem=''.join(nlist)
            #nitem=item.split('|')[0].split('.')[0]
            nitem=item
            if nitem in dict1:
                iseq=dict1[nitem]
                out1.write('>%s\n%s'%(nitem,iseq))
                i+=1; f+=1
        if i>0:
            out1.write('##\n')
            x+=1
    line2=file2.readline()
file2.close()
out1.close()
print "No. of items written: ", f
print "No. of chunks written: ", x
print "Done!"
            
            
    
