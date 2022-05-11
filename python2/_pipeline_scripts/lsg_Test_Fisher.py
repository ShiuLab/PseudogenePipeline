"""
python Test_Fisher.py file 0 #without q-value
  
python Test_Fisher.py file 1 # with q-value
05 08 2008, Revice it to two tail

0800805 add a choice about calculate q-value or not



"""

#To run this piece of garbage the 1st parameter is the input file 5col file in
#the format:
#category  group1+  group1-  group2+  group2-       where - means does not have
#this category and + means does have this category.  The 2nd parameter is 1
#Additionally load module R before running and make sure there are no special
#chars in the input file like ' or ,

print "USE THIS SCRIPT IF FISHERINP HAS 5 COLUMNS (NAME-v1-v2-v3-v4)"
print "If there's an error, change the following line in the script"
print "Line 76: if len(L)==5: "
def fisher(tuple_t):# every 4 number tuple

	"""
    * k number of objects in the selection having a given property
    * n size of the selection
    * C number of objects in the population having this property
    * G size of the population
	"""
	t1=int(tuple_t[0])
	t2=int(tuple_t[1])
	t3=int(tuple_t[2])
	t4=int(tuple_t[3])

	f = FisherExactTest()
	p=f.pvalue(t1,t2+t1,t1+t3,t1+t2+t4+t3)
	en=f.enrichment(t1,t2+t1,t1+t3,t1+t2+t4+t3)
	
	if en >= 1:
		en="+" #using right hand side Ha :u
		
	else:
		en="-"
		
	#print p
	
	pv=p[2]
	#print pv
	return en,pv

def get_q(file,path):
	print "Here"
	print file
	tmpR=open("%s_tmp.R" % file,"w")
	tmpR.write("setwd(%r)\n" % ("%s"  % path))
	tmpR.write("source(%r)\n" % "/home/hufnag30/bin/qvalue/R/qvalue.R")
	tmpR.write("test<-read.table(%r,header=F)\n" % ("%s" % file))
	tmpR.write("Q<-qvalue(test$V7,pi0.method=%r,lambda=0)\n" % "bootstrap")
	tmpR.write("q<-Q$qvalue\n")
	tmpR.write("test$V8<-q\n")
	tmpR.write("write.table(test,file=%r,quote=F,eol= %r ,row.names=F,col.names=F,sep=%r)\n" % ('%s.qvalue' % file,'\n','\t'))
	tmpR.close()
	#sys.exit()
	os.system("R CMD BATCH %s_tmp.R" % file)
	#print "Done"
	#sys.exit()


if __name__ == '__main__':
	
	import os,sys
	#sys.path.append("/home/zou/pycodes")
	from fisher import FisherExactTest
	file=sys.argv[1]
	q=sys.argv[2] # 1 or 0 1: calculate the q-value,0:do not calculate
        genecount=0
	path = os.getcwd()
        
	temp_fisher=open("%s_temp_fisher.test" % file,"w")
	#new_temp_fisher=open("%s.fisher.pqvalue.withgenenames" % file,"w")
        nlist=[]
	for line in open(file,"r"):
                if line.startswith('#')==True:
                        #line=line.strip()
                        #temp_fisher.write('%s\tresult[0]\tresult[1]\n')
                        pass
                else:
                        L=line.strip().split("\t")
        		if len(L)==5:                       
        
        			"""
        			use_chi=1
        			for i in L[1:]:
        				#print i
        				#print type(i)
        				if int(i) < 5:
        					use_chi=0
        					break
        			if use_chi==1:
        				temp_chi.write("%s" % line)
        			else:
        				result = fisher(L[1:])
        			"""
        			result = fisher(L[1:])
                        
        			#new_temp_fisher.write("%s\t%s\t%s\t%s\n" % \
                                    #              (L[0],"\t".join(L[1:]),result[0],result[1]))
        			temp_fisher.write("%s\t%s\t%s\t%s\n" % \
                                                  (L[0],"\t".join(L[1:]),result[0],result[1]))
        			
                        if L[0] not in nlist:
                                #print L[0]
                                genecount+=1
                                nlist.append(L[0])
	temp_fisher.close()
	#print "File closed"
	#sys.exit()
	#new_temp_fisher.close()
	print "No. of items analysed: ", genecount
	if q=="0":
		os.system("mv %s_temp_fisher.test %s.fisher.pqvalue" % (file,file))
		
	if q=="1":
		#add the chi result to fisher result
		#os.system("perl /home/hanada/bin/Chi_test.pl temp_chi>temp")
		#os.system("python /home/shiu/codes/FileUtility.py -f -c 3 -i temp_chi_test")
		get_q("%s_temp_fisher.test" % file,path)
		os.system("mv %s_temp_fisher.test.qvalue %s.fisher.pqvalue" % (file,file))

		#os.system("cat temp_fisher.test.qvalue temp_chi_test > %s.Chisq_fisher.pqvalue" % file)
		#os.system("mv temp_fisher.test.qvalue %s.fisher.pqvalue" % file )
		#os.system("mv temp_chi_test %s.chisq.pqvalue" % file)
                #sys.exit()
		os.system("rm %s_temp_fisher.test" % file)
		os.system("rm %s_temp_fisher.test_tmp.R" % file)
		os.system("rm %s_temp_fisher.test_tmp.Rout" % file)
		#os.system("rm temp_fisher.test")

print "Now filtering the output file..."
file2=open(file+".fisher.pqvalue",'r')
out1=open(file+".fisher.pqvalue.sign+pq005.tab",'w')
out2=open(file+".fisher.pqvalue.sign+pq005.tab.morethan10genes.tab",'w')
out1.write('#python %s\n'%(' '.join(sys.argv)))
out2.write('#python %s\n'%(' '.join(sys.argv)))
line2=file2.readline()
while line2:
        if line2.startswith('#'):
                out1.write(line2)
        else:
                tab2=line2.strip().split('\t')
                #print tab2
                if tab2[5]=='+' and float(tab2[7])<0.05:
                        print "in"
                        out1.write(line2)
                        if int(tab2[1])>=10:
                                out2.write(line2)
        line2=file2.readline()
file2.close()
out2.close()
                
        
        
