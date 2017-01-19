import sys
print "INP1: INP file"
print "INP2: table/fasta (Doesnt work for fasta)"
print "INP3: Suffix to be used for marking (lyr/thal...)"
file1=open(sys.argv[1],'r')
line1=file1.readline()
ftp=sys.argv[2].strip()
out1=open(sys.argv[1]+".repmarked",'w')
dict1={}
list1=[]
suf=sys.argv[3]
m=0
if ftp=='table':
    while line1:
        if line1.startswith('#')==True:
            if m==0:
                pass
            else:
                #print list1
                if len(list1)==0:
                    pass
                else:
                    for line in list1:
                        tab1=line.strip().split('\t')
                        g1=tab1[0]
                    if g1 not in dict1:
                        dict1[g1]=1
                        for line in list1:
                            #out1.write(line)
                            tab1=line.strip().split('\t')
                            g1=tab1[0]
                            count=1
                            out1.write('%s|%s%s\t%s\n'%(g1,suf,count,'\t'.join(tab1[1:])))
                    else:
                        count=dict1[g1]
                        count+=1
                        dict1[g1]=count
                        for line in list1:
                            tab1=line.strip().split('\t')
                            g1=tab1[0]
                            out1.write('%s|%s%s\t%s\n'%(g1,suf,count,'\t'.join(tab1[1:])))
                    list1=[]
                    out1.write('##\n')
        else:
            list1.append(line1)
        m+=1
        line1=file1.readline()
        
       
elif ftp=='fasta':
    sys.exit()
    while line1:
        if line1.startswith('>')==True:
            name=line1[1:-1]
            if name not in dict1:
                dict1[name]=1
                out1.write(line1)
            else:
                count=dict1[name]
                count+=1
                dict1[name]=count
                out1.write('>%s-%s\n'%(name,count))
        else:
            out1.write(line1)
        line1=file1.readline()
file1.close()
out1.close()


