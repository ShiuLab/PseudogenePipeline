#This script is desinged to run the Merge EST pipeline in a wrapper provided
#that the proper modules are loaded ahead of time (Vmatch, EMBOSS, RepeatMasker).
#Created by David E. Hufnagel on Aug 27, 2012
#WARNING: RUN AS A JOB IT SHOULD TAKE ~1 DAY

import sys, os

base = sys.argv[1]  #the name of the input ESTs file
crap = sys.argv[2]  #the file containing contaminats DNA for extraction
proc = sys.argv[3]  #the number of processers to use

#os.system("python ~/Shiu/Scripts/FastaManagerDave.py -f del_redun_by_seq -fasta %s" % (base))
#base += ".mod"
os.system("mkvtree -db %s -dna -pl -allout -v" % (crap))
os.system("vmatch -qmaskmatch X -d -p -l 50 -exdrop 1 -identity 90 -q %s %s  > output.txt" % (base, crap))
os.system("trimest %s %s.noPolA" % (base, base))
os.system("mkvtree -db %s.noPolA -dna -pl -allout -v" % (base))
os.system("vmatch -d -p -l 50 -exdrop 1 -identity 95 -v -showdesc 0 -dbcluster 100 10 -nonredundant %s.noPolA.nonRed %s.noPolA > output.info" % (base, base))
os.system("python ~/Shiu/Scripts/NameChangeManager.py -f FastaKillSpace2 -inp %s.noPolA.nonRed -out %s.noPolA.nonRed.mod" % (base, base))
os.system("python ~/Shiu/Scripts/FastaManagerDave.py -f size_filter -fasta %s.noPolA.nonRed.mod -T 100" % (base))
##os.system("RepeatMasker -cutoff 200 -par %s -xsmall -species viridiplantae -gff %s.noPolA.nonRed.mod_T100.fa" % (proc, base))
##os.system("python ~/Shiu/Scripts/parse_RepeatMasker_gff.py %s.noPolA.nonRed.mod_T100.fa.out 225 30" % (base))
##os.system("python ~/Shiu/Scripts/1_6_4_4_removeSeqsFromFasta.py %s.noPolA.nonRed.mod_T100.fa %s.noPolA.nonRed.mod_T100.fa.out.225Cutoff30.0divg.4col %s.noPolA.nonRed.mod_T100.fa.RMfilt" % (base, base, base))
##os.system("python ~/Shiu/Scripts/Seq_specialCharToNnt.py %s.noPolA.nonRed.mod_T100.fa.RMfilt %s.noPolA.nonRed.mod_T100.fa.RMfilt.mod" % (base, base))
##os.system("mkdir 1_FilterOutCrap")
##os.system("mkdir 2_TrimPolyAs")
##os.system("mkdir 3_RemoveDupsWithVmatch")
##os.system("mkdir 4_RepeatMasking")
##os.system("mv output.txt %s 1_FilterOutCrap" % (base))
##os.system("mv %s.noPolA* 2_TrimPolyAs" % (base))
##os.system("mv 2_TrimPolyAs/%s.noPolA.nonRed* ." % (base))
##os.system("mv %s.noPolA.nonRed* output.info 3_RemoveDupsWithVmatch" % (base))
##os.system("mv 3_RemoveDupsWithVmatch/%s.noPolA.nonRed.mod_T100.fa.* ." % (base))
##os.system("mv %s* 4_RepeatMasking" % (base))
##print "\nNow change the name so it isn't so long (it's the ...4_RepeatMasking/... .RMfilt file)\n"
