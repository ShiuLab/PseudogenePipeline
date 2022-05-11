#This script is designed run all the scripts involved in parsing the GMAP output
#into it's final format
#Created by David E. Hufnagel on Aug 29, 2012
#WARNING: MANY PARAMETERS AREN'T VARIABLES IN THIS WRAPPER

import sys, os

gmap = sys.argv[1]     #the name for the gmap output
database = sys.argv[2] #the input database for GMAP
idt = sys.argv[3]      #the %ID for GMAP parsing
cov = sys.argv[4]      #the coverage for GMAP parsing

os.system("python ~/Shiu/Scripts/gmap_chunks_step1.py %s" % (gmap))
os.system("python ~/Shiu/Scripts/gmap_chunks_step1.5_IDTCOVparse.py %s.breaks %s %s" % (gmap, cov, idt))
os.system("python ~/Shiu/Scripts/gmap_chunks_step3.py %s.breaks.%sidt%scovfilt" % (gmap, idt, cov))
os.system("python ~/Shiu/Scripts/mark_repeatedID_col0_chunks.py %s.breaks.%sidt%scovfilt.4col table thal" % (gmap, idt, cov))
os.system("python ~/Shiu/Scripts/chunks_get_minmax.py %s.breaks.%sidt%scovfilt.4col.repmarked" % (gmap, idt, cov))
os.system("grep -v mitochondria %s.breaks.%sidt%scovfilt.4col.repmarked | grep -v chloroplast > %s.breaks.%sidt%scovfilt.4col.repmarked.noCM" % (gmap, idt, cov, gmap, idt, cov))
os.system("grep -v mitochondria %s.breaks.%sidt%scovfilt.4col.repmarked.minmax | grep -v chloroplast > %s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM" % (gmap, idt, cov, gmap, idt, cov))
os.system("python ~/Shiu/Scripts/chunks_single_multiple_getstretch4_fast_type1.py %s.breaks.%sidt%scovfilt.4col.repmarked.noCM %s" % (gmap, idt, cov, database))
os.system("python ~/Shiu/Scripts/get_longest_isoform.py %s.breaks.%sidt%scovfilt.4col.repmarked.noCM.fa thal" % (gmap, idt, cov))
os.system("python ~/Shiu/Scripts/1_6_5b_getLongestFromMinmax.py %s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM %s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM.longest" % (gmap, idt, cov, gmap, idt, cov))
os.system("python ~/Shiu/Scripts/1_6_6_longestToBlast.py %s.breaks %s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM.longest %s.breaks.%sidt%scovfilt.4col.repmarked.minmax.noCM.longest.blast" % (gmap, gmap, idt, cov, gmap, idt, cov))
print "\ndone with wrapper!\n"
