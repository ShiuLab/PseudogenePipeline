#This script is designed to use Shinhan's scripts to obtain retroPseudogenes
#given that pseudogenes are already acquired

import sys, os

gff = sys.argv[1]        #genomic gff file
pseu = sys.argv[2]       #fasta file of pseudogenes
prot = sys.argv[3]       #Translated cds fasta file
pseu4col = sys.argv[4]   #4col file of pseudogenes
scriptDir = sys.argv[5]  #the directory containing the scripts this wrapper calls
scriptDir = "/" + scriptDir.strip("/") + "/"

os.system("python %s5_2_5_getLongestCDSgff.py %s %s.CDS" % (scriptDir, gff, gff))
os.system("python %sTabularManager.py -f shiftCols -cols 1,4,5,7,8,9 -inp %s.CDS -out %s.cds.coord" % (scriptDir, gff, gff))
os.system("python %sscript_6.4.3_prot_ejunctionDave.py %s.cds.coord" % (scriptDir, gff))
os.system("cat %s %s > pseuProtTemp.aa.fa" % (pseu, prot))
os.system("python %s1_3_E_make2colHomologyOfPseus.py %s %s.PSvFG.2col" % (scriptDir, pseu4col, pseu4col))
os.system("python %sBlastUtilityDave.py -f batch_bl2 -i pseuProtTemp.aa.fa -p blastx -g %s.PSvFG.2col -d 1 -W 3 -bdir ~shius/bin/blast/bin" % (scriptDir, pseu4col))
os.system("python %sParseBlast.py -f parse_align2 -blast %s.PSvFG.2col.bl2seq.raw" % (scriptDir, pseu4col))
os.system("python %sParseBlast.py -f parse_gap -blast %s.PSvFG.2col.bl2seq.raw" % (scriptDir, pseu4col))
#os.system("python %s5_2_4_replacePACIDinGap.py %s.cds.coord.exon_junc %s %s.cds.coord.exon_junc.rlNames" % (scriptDir, gff, gff, gff))
#os.system("python %sscript_6.4.3_retro_or_not.py %s.cds.coord.exon_junc.rlNames %s.PSvFG.2col.bl2seq.raw.gap %s.PSvFG.2col.bl2seq.raw.mod > %s.PSvFG.2col.bl2seq.raw.log" % (scriptDir, gff, pseu4col, pseu4col, pseu4col))
#os.system("" % ()) #make the actual list of retropseudogenes

#Remove not needed files
os.system("mkdir Temp")
os.system("mv temp Temp")
os.system("mv pseuProtTemp.aa.fa Temp")
os.system("mv %s.CDS Temp" % (gff))
os.system("mv %s.cds.coord.exon_junc Temp" % (gff))
