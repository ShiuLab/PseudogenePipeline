#This script is designed run all the scripts involved in parsing the GMAP output
#into it's final format
#Created by David E. Hufnagel on Aug 29, 2012
#WARNING: MANY PARAMETERS AREN'T VARIABLES IN THIS WRAPPER

import sys, os
print "INP1: GMAP output file"
print "INP2: FASTA of database used for GMAP"
print "INP3: IDT threshold"
print "INP4: COV threshold"
print "INP5: Suffix you want to attach to the output

gmap = sys.argv[1]     #the name for the gmap output
db = sys.argv[2] #the input database for GMAP
idt=sys.argv[3]; cov=sys.argv[4]; suf=sys.argv[5]
dir1='/home/moghegau/scripts/projects/3_RNA_genes/GMAP_scripts'

os.system("python %s/gmap_chunks_step1.py %s" % (dir1,gmap))
os.system("python %s/gmap_chunks_step1.5_IDTCOVparse.py %s.breaks %s %s" \
          % (gmap,idt,cov))
os.system("python %s/gmap_chunks_step3.py %s.breaks.%sidt%scovfilt" % \
          (gmap,idt,cov))
os.system("python %s/mark_repeatedID_col0_chunks.py " \
          "%s.breaks.%sidt%scovfilt.4col table %s" % (gmap,idt,cov,suf))
os.system("python %s/chunks_get_minmax.py " \
          "%s.breaks.%sidt%scovfilt.4col.repmarked" % (gmap,idt,cov))
os.system("grep -v mitochondria %s.breaks.%sidt%scovfilt.4col.repmarked " \
          "| grep -v chloroplast > %s.breaks.%sidt%scovfilt.4col.repmarked.noCM" \
          % (gmap, idt, cov, gmap, idt, cov))
os.system("grep -v mitochondria %s.breaks.%sidt%scovfilt.4col.repmarked.minmax " \
          "| grep -v chloroplast > %s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM" \
          % (gmap, idt, cov, gmap, idt, cov))
os.system("python %s/chunks_single_multiple_getstretch4_fast_type1.py " \
          "%s.breaks.%sidt%scovfilt.4col.repmarked.noCM %s" % \
          (gmap, idt, cov, database))
os.system("python ~/Shiu/Scripts/get_longest_isoform.py " \
          "%s.breaks.%sidt%scovfilt.4col.repmarked.noCM.fa %s" % \
          (gmap, idt, cov, suf))
print "\nDone with wrapper!\n"
