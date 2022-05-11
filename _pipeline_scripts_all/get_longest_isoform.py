import sys
'''
From lyrata GMAP sequence file, get the lyrortholog with longest sequence
'''
print "INP1: Repmarked FASTA file"
print "INP2: Suffix used in the repmarked FASTA file (eg: lyr/thal)"
file1=open(sys.argv[1],'r').readlines()
out1=open(sys.argv[1]+".longest",'w')
out1.write('#python %s\n'%(' '.join(sys.argv)))
suf=sys.argv[2]

dict1={}
for i in range(0,len(file1)):
    line=file1[i]
    if line.startswith('>'):
        name=line[1:].strip().split('|')[0]
        seq=file1[i+1].strip()
        if name not in dict1:
            dict1[name]=seq
        else:
            if len(seq)>len(dict1[name]):
                dict1[name]=seq

for name in dict1:
    seq=dict1[name]
    out1.write('>%s|%s\n%s\n'%(name,suf,seq))
out1.close()    
