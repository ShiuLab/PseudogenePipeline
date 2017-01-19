import sys, readFile
'''
AthAly000       Chr1    930828  3255693
AthAly000       scaffold_1      1125712 3818257
AthAly001       Chr1    4968142 6444076
AthAly001       scaffold_1      6134005 7852486
AthAly002       Chr1    9603178 10647511
AthAly002       scaffold_1      12120807        13428671
'''

print "INP1: All features 4col file"
print "INP2: Orthopairs 2col file"
print "INP3: What name you want to give to the pair (index0)"

file1=open(sys.argv[1],'r')
line1=file1.readline()
name=sys.argv[3]
dict1={}
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        g1=tab1[0]; rest='\t'.join(tab1[1:])
        if g1 not in dict1:
            dict1[g1]=rest
        else:
            print "Gene repeat: ", g1
    line1=file1.readline()
file1.close()

file2=open(sys.argv[2],'r')
out1=open(sys.argv[2]+".links",'w')
line2=file2.readline()
link=0
while line2:
    if line2.startswith('#'):
        pass
    else:
        tab2=line2.strip().split('\t')
        g1=tab2[0]; g2=tab2[1]
        if g1 in dict1 and g2 in dict1:
            link+=1
            r1=dict1[g1]; r2=dict1[g2]
            out1.write('%s%s\t%s\n%s%s\t%s\n'%(name,link,r1,name,link,r2))
    line2=file2.readline()
file2.close(); out1.close()
print "Links written to OUT: ", link
print "Done!"
            
        
            
                                 
