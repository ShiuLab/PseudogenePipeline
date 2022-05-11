#This script is designed to get pseudogene CDS's from standard pseudogene files
#Created by David E. Hufnagel on Jan 14, 2013
#Updated on Mar 19, 2013 to handle proper .ref files

import sys, os

ref = sys.argv[1]       #input reference file connecting original names to big names (alternative names will be developed from the big names)
dis = sys.argv[2]       #input disable_count pseudogene file
four = sys.argv[3]      #input 4col pseudogene file with code names
genome = sys.argv[4]    #input genome file
scriptDir = sys.argv[5] #the directory where the scripts are located

scriptDir = "/" + scriptDir.strip("/") #get rid of '/' after scriptDir


os.system("python %s/FastaManagerGaurav.py -f get_stretch4 -fasta %s -coords %s" % (scriptDir, genome, four))
fasta = "%s.fa" % (four)
os.system("python %s/1_3_1_C_makeAltPseuName.py %s %s %s.alt" % (scriptDir, ref, fasta, fasta))
os.system("python %s/script_step7Dave.py %s.alt %s > step7.log" % (scriptDir, fasta, dis))
os.system("python %s/script_step7.1.py %s.ps_cds " % (scriptDir, dis))
os.system("python %s/1_3_1_C_undoAltPseuName.py %s %s.ps_cds.rid_allstop %s.ps_cds.rid_allstop.origName" % (scriptDir, ref, dis, dis))
os.system("mkdir Temp")
#os.system("mv %s.ps_cds* Temp" % (dis))
#os.system("mv Temp/%s.ps_cds.rid_allstop.origName ." % (dis))
os.system("mv %s.alt Temp" % (fasta))
os.system("mv %s.ps_cds.rid_allstop.origName %s.CDS" % (dis, fasta))
