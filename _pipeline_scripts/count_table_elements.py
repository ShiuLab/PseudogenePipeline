#how many different elements r present in a table column u specify

import sys
print "Assumption is this is a tab-delimited file"
print "INP1: Table name"
print "INP2: Index number of column to scan"

#print "Output: *.idlist"
#print "Column name.startswith(AT)==True"
file1=open(sys.argv[1], 'r') #table file name
col=int(sys.argv[2]) #the INDEX, not column
#out1=open(sys.argv[1]+".idlist", 'w')
line1=file1.readline()
list1=[]
dict1={}
m=0
s=0
repeat=0
while line1:
    if line1.startswith('#')==True:
        pass
    else:
        tab1=line1.strip().split('\t')
        item=tab1[col]        
        #item21=tab1[col].split('-')
        #item2=item21[0]
        #item1=item2.split('.')
        #item=item21[0]
        #item=item2[0]
        if item not in dict1:# and item.startswith('>PUT'):            
            #print item
            dict1[item]=1
            #out1.write('%s\t%s\n'%(item,tab1[1]))
            m+=1
            #if item.startswith('AT')==True:
             #   dict1[item]=1
              #  m+=1
        else:
            repeat+=1
    s+=1
    if s%1000000==0:
        print s
    line1=file1.readline()
print "Number of unique reads: ", len(dict1.keys())
print "Number of repeated reads: ", repeat
print "No. of different items in INDEX-",col,"-are: ", m
