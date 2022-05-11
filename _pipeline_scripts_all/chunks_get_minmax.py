from __future__ import division
import sys
print "INP1: Chunks GFF file"
file1=open(sys.argv[1], 'r')
line1=file1.readline()
dict1={}
m=0
n=0
d=0

out1=open(sys.argv[1]+".minmax",'w')
list1=[]
while line1:
    if line1.startswith('##')==True:
        if len(list1)==1:
            for line in list1:
                out1.write(line)
        else:
            list2=[]
            for line in list1:
                tab1=line.strip().split('\t')
                g1=tab1[0]
                chr1=tab1[1]
                st=int(tab1[2])
                end=int(tab1[3])
                list2.append(st)
                list2.append(end)
            
            nst=min(list2)
            nend=max(list2)
            out1.write('%s\t%s\t%s\t%s\n'%(g1,chr1,nst,nend))
        list1=[]
        #out1.write('##\n')
    elif line1.startswith('#'):
        pass
    else:
        list1.append(line1)
    line1=file1.readline()
