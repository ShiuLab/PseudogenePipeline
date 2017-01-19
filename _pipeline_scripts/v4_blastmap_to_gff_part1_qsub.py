#split table and run blastmap
from __future__ import division
import sys
import os
print "#######"
print "INP1:   The TEST 4col File"
print "INP2:   How many divisions wanted?"
print "INP3:   Reference 4col file"
print "#######"
file1=open(sys.argv[1], 'r') #name of PUT GFF file
fid=sys.argv[1]
pwd=os.getcwd()
by=int(sys.argv[2])  #how many pieces
pwd=os.getcwd()  #P. Working directory PWD
db=sys.argv[3]  #Database file name (GFF of all genes)

line1=file1.readline()
dict1={}
m=1
print "Here"
os.system('python /home/moghegau/scripts/Fastaformat2.py -f split_table -table %s -n %s'%\
          (fid,by))
out2=open('jobv4_commands.txt','w')
i=1
while i<(by+2): #$$$$$ change to by+1
    out2.write('python /home/hufnag30/Shiu/Scripts/v4_blastmap_to_gff_part2_smallchange.py ' +\
               '%s %s.frag%s n\n' % (db,fid,i))
    '''
    out2.write('python ~moghe/lab_people/melissa/sorf2/' +\
               'v4_blastmap_to_gff_part2_smallchange.py ' +\
               '%s %s.frag%s y 25\n' % (db,fid,i))
    '''    
    i+=1
out2.close()
print "Done!"
