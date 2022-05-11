import sys
file1=open(sys.argv[1], 'r')
out1=open(sys.argv[1]+".fa", 'w')
m=0
n=0
line1=file1.readline()
def rmlb(line):
    if line.endswith('\n'):
        line=line[:-1]
    else:
        pass
    return line

while line1:
    if line1.startswith('#python'):
        pass
    
    elif line1.startswith('>')==True:        
        if m==0:
            pass
        else:
            #print line1
            seq1=''.join(seq)
            if len(seq1)==0:                
                n+=1                
            else:
                out1.write('>%s\n%s\n'%(gname, seq1))            
        gname=line1[1:-1]; seq=[]
        
    elif line1.startswith('#')==True:            
            seq1=''.join(seq)
            if len(seq1)==0:                
                n+=1                
            else:
                out1.write('>%s\n%s\n'%(gname, seq1))
            gname=''; seq=[]
            out1.write(line1)
            m=0            
    elif line1.startswith('==='):
        pass    
    else:   
        line2=line1.strip()
        if line2!='':
            seq.append(line2)
    
    m+=1
    if m%1000000==0:
	print m
    line1=file1.readline()
        
seq1=''.join(seq)
#print seq1
if len(seq1)==0:
    n+=1
    print "Here and wrong"
else:
    out1.write('>%s\n%s\n'%(gname, seq1))
    print('%s\t%s\n'%(gname, len(seq1)))
file1.close()
print "Number of zero-length sequences, now skipped: ", n
    
