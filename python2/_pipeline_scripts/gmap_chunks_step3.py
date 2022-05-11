import sys
print "INP1: Output of step2"
file1=open(sys.argv[1],'r')
line1=file1.readline()
dict1={}
m=0
n=0
out1=open(sys.argv[1]+".4col",'w')
#out2=open(sys.argv[1]+".4col.mod",'w')
#use single_multiple_getstretch4.py
while line1:
    if line1.startswith('#')==True:
        out1.write(line1)
    else:
        tab1=line1.strip().split('\t')        
        g1=tab1[0]
        '''
        g1sp=g1.split('|')
        csp=g1sp[0].split('_')
        chr1=csp[1]
        lsp=g1sp[3].split('-')
        sign=g1sp[1]
        if sign=='+':
            cl=lsp[0]
            cr=lsp[1]
        elif sign=='-':
            cl=lsp[1]
            cr=lsp[0]
        #cl=lsp[0]
        #cr=lsp[1]
        
        '''
        #This is actually the right loop. Am altering it for sORF sequences"
        loc=tab1[8].split(':')
        chr1=loc[0]
        span=loc[1].split('-')
        cl=span[0]
        cr=span[1]
        
        out1.write('%s\t%s\t%s\t%s\n'%(g1,chr1,cl,cr))
        #out2.write('%s\t%s\t%s\t%s\n'%(g1,chr1,cl,cr))
    line1=file1.readline()
file1.close()
out1.close()
print "Done!"
 
