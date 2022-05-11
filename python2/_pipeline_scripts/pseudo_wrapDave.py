###############################################################################
#
# 02/16/12 Shinhan - Bugs... When file is empty, there is no warning. Also,
#          fast3 version change so the program name is changed.
#
# 11/21/11 Shinhan - Modify the pipeline script so it can better be used for
#          standalone applications. Make p_codes, o_codes, f_dir, and blosum
#          as required parameters. Also, activate the blast filtering part of
#          the codes.
#
# 10/14/10 Shinhan - Pseudogene pipeline wrapper. Assuming the BLAST run is
#          finished. Require a parameter file to run.
# 
#

import sys, os

if len(sys.argv) != 2:
	print "Usage: pseudo_wrap.py parameter_file"
	print " Parameter file includes the following required values:"
	print "   b_out     Blast output (protein vs. genome)"
	print "   p_seq     Query protein sequence file"
	print "   g_seq     Genome/DNA sequence"
	print "   il_t      Intron length threshold"
	print "   b_filter  Filter blast output (y) or not (n)"
	print "   p_codes   Path to pseudo pipeline scripts"
	print "   o_codes   Path to blast parser script"
	print "   f_prog    Fasta 3.x executable for tfasty"
	print "   blosum    Path to BLOSUM matrix"
	print " And the following optional values"
	print "   ev_t      Evalue threshold in -log(E) format, default 5"
	print "   id_t      Identity threshold in %, default 40"
	print "   ml_t      Matching region length threshold in aa, default 30"
	print "   ml_p      Matching region proprotion to query threshold, default 0.05"
	print "   "
	print " Parameter file example - /home/public/script/pseudo/example/parameter"
	sys.exit(0)

###
# Read and check parameters
###
print "Check parameter file..."
inp = open(sys.argv[1])
inl = inp.readlines()
P   = {}
for i in inl:
	L = i.strip().split("=")
	if len(L) != 2:
		print "Parameter file format problem:",[L]
		print "Quit!"
		sys.exit(0)

	P[L[0]] = L[1]


# Check if required values are missing
required = ["b_out","p_seq","g_seq","il_t","p_codes",
						"b_filter","o_codes","f_dir","f_prog","blosum"]
p_missing = 0
for i in required:
	if i not in P:
		print " missing:",i
		p_missing = 1

if p_missing:
	sys.exit(0)
	
# Check if blast output and sequence files exist
f_missing = 0
r_files = [P["b_out"],P["p_seq"],P["g_seq"]]
for i in r_files:
	if not os.path.isfile(i):
		print " file does not exist:",i
		f_missing = 1

if f_missing:
	sys.exit(0)
	
p_codes = P["p_codes"]
o_codes = P["o_codes"]
f_dir   = P["f_dir"]
f_prog  = P["f_prog"]
blosum  = P["blosum"]

# Check if other values are set, if not default is used.
# Set default values
default = {"ev_t"    :"5",
					 "id_t"    :"40",
					 "ml_t"    :"30",
					 "ml_p"    :"0.05"}
for i in default:
		if i not in P:
			P[i] = default[i]
			print " default: %s=%s" % (i,P[i])

###
# Filter BLAST output
###
f_flag = 0
if "b_filter" in P and P["b_filter"] == "y":
	f_flag = 1
	try:
		print "Filter %s..." % P["b_out"]
		print " E:%s I:%s L:%s P:%s" % (P["ev_t"],P["id_t"] ,P["ml_t"] ,P["ml_p"])
	except KeyError:
		print "Missing filtering parameter(s)! Quit!"
		sys.exit(0)
	os.system("python %s/ParseBlast.py -f get_qualified4 -blast " % o_codes +\
					  "%s -fasta %s -E %s -I " % (P["b_out"],P["p_seq"],P["ev_t"])  +\
					  "%s -L %s -P %s -Q 1 "   % (P["id_t"] ,P["ml_t"] ,P["ml_p"])  +\
					  "> log_step1")
					  
	ml_p = int(float(P["ml_p"])*100)
	os.system("mv %s_E%sI%sL%sP%sQ1.qlines %s_parsed" % \
						(P["b_out"],P["ev_t"],P["id_t"],P["ml_t"],ml_p,P["b_out"]))
	P["b_out"] = "%s_parsed" % P["b_out"]
	
	# Check if there is problem with the file


###
# Get pseudoexon and phase 1 pseudogene
###
print "Get pseudoexons..."
os.system("python %s/script_step2e.py %s %s > log_step2_pe" % \
																								(p_codes,P["b_out"],P["il_t"]))
pe = "%s_G%s.PE" % (P["b_out"],P["il_t"])
print " pseudoexon file:",pe

print "Get phase 1 pseudogene..."
os.system("python %s/script_step3bDave.py %s %s > log_step3_phase1_ps" % \
																								(p_codes,pe,P["il_t"]))
ps1 = "%s_I%s.PS1" % (pe,P["il_t"])
print " phase1 ps file:",ps1

#
# Prepare ps1 sequnce file and pair file
###
print "Get pair file and subject coordinates..."
os.system("python %s/script_step3.5.py %s" % (p_codes,ps1))
print " pair file  : %s.pairs" % ps1
print " coordinates: %s.subj_coord" % ps1

# Get ps1 sequences
print "Get phase 1 pseudogene sequences..."
ps1coord = ps1 + ".subj_coord"
os.system("python %s/FastaManager.py -f get_stretch4 -fasta " % o_codes + \
		  "%s -coords %s > log_step4_phase1_ps_seq" % (P["g_seq"],ps1coord))
print " phase 1 ps sequence: %s.subj_coord.fa" % ps1

###
# Run SW
###
print "Find stop and framshifts..."
os.system("python %s/BlastUtility.py -f batch_sw -p %s -g " % (o_codes,f_prog)+\
		  "%s.pairs -i %s -j %s.fa -bdir " % (ps1,P["p_seq"],ps1coord)    + \
		  "%s -d 1"         % f_dir)
print " Smith-Waterman outputs: %s_pairs.sw.*" % ps1
	
os.system("python %s/script_step6.py %s %s.pairs.sw.out" % \
																				(p_codes,blosum,ps1))
print" Final output: %s_pairs.sw.out.disable_count" % ps1

print "\nThe pseudogene pipeline has finished!"







