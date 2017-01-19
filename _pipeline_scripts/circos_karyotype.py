import sys, readFile
'''
Make Circos karyotype file
#Arabidopsis Thaliana drawn in blue
chr     -       At1     At1     0       30427670        blue
chr     -       At2     At2     0       19698288        blue
chr     -       At3     At3     0       23459829        blue
chr     -       At4     At4     0       18585055        blue
chr     -       At5     At5     0       26975501        blue
band    At1     AT1G01010       AT1G01010       3631    5899    blue
band    At1     AT1G01020       AT1G01020       5928    8737    blue
band    At1     AT1G01030       AT1G01030       11649   13714   blue
band    At1     AT1G01040       AT1G01040       23146   31227   blue
band    At1     AT1G01046       AT1G01046       28500   28706   blue
'''
print "INP1: Genome FASTA file"
print "INP2: 4col file"
print "INP3: Color you want to give to this genome (blue/red/green/black)"
print "INP4: Organism name"

dict1={}; dict2={}
dict1=readFile.readFasta(sys.argv[1],dict1)
col=sys.argv[3]; orgName=sys.argv[4]

out1=open(sys.argv[1]+".karyotype",'w')
out1.write('#%s drawn in %s\n'%(orgName,col))
for name in dict1:
    if 'mito' in name or 'chloro' in name or 'Scaffold' in name or 'ChrM' in name:
        pass
    else:
        seqSize=len(dict1[name])
        out1.write('chr\t-\t%s\t%s\t0\t%s\t%s\n'%(name,name,seqSize,col))
        dict2[name]=1        
##        try:
##            if int(name)<=8:
##                seqSize=len(dict1[name])
##                out1.write('chr\t-\t%s\t%s\t0\t%s\t%s\n'%(name,name,seqSize,col))
##                dict2[name]=1
##        except:
##            seqSize=len(dict1[name])
##            out1.write('chr\t-\t%s\t%s\t0\t%s\t%s\n'%(name,name,seqSize,col))
        
                
    
file1=open(sys.argv[2],'r') #4-col file
line1=file1.readline()
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        g1=tab1[0]; chr1=tab1[1]; st=tab1[2]; end=tab1[3]
        if chr1 in dict2:
            out1.write('band\t%s\t%s\t%s\t%s\t%s\t%s\n'%\
                       (chr1,g1,g1,st,end,col))
        else:
            print "Chrom. not present in FASTA file: ", chr1
            #sys.exit()            
    line1=file1.readline()
file1.close(); out1.close()

print "Done!"


    
