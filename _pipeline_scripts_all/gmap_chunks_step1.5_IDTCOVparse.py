import sys
print "INP1: Chunks Breaks file. Output of gmap_chunks_part1.py"
print "INP2: Coverage"
print "INP3: Identity"
file1=open(sys.argv[1], 'r')

line1=file1.readline()
cov=float(sys.argv[2].strip())
idt=float(sys.argv[3].strip())
out1=open(sys.argv[1]+".%sidt%scovfilt"%(sys.argv[3].strip(),sys.argv[2].strip()), 'w')
m=0
true=0
g=0
while line1:
    if line1.startswith('#')==True:
        if true==0:
            out1.write(line1)
            true=1
    else:
        tab1=line1.strip().split('\t')
        pname=tab1[0]
        nhits=tab1[1]
        cov1=float(tab1[3])
        idt1=float(tab1[4])        
        if cov1>=cov and idt1>=idt:
            #print pname,cov1,idt1
            out1.write(line1)
            m+=1
            true=0
        g+=1
        #sys.exit()
    line1=file1.readline()
print "Total number of lines: ", g
print "Written to output: ", m
        
